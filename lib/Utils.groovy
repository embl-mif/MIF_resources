import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.security.MessageDigest

class Utils {
    static String serializeMobieInputUri(def image, Path tableDir) {
        def inputUri
        if (image instanceof java.nio.file.Path) {
            def imageAbs = image.toAbsolutePath().normalize()
            try {
                inputUri = tableDir.relativize(imageAbs).toString()
            } catch (Exception _ignored) {
                // Fall back to absolute URI when paths are not relativizable (e.g. different roots).
                inputUri = imageAbs.toString()
            }
        } else {
            // URI inputs are already the desired serialized representation.
            inputUri = image.toString()
        }

        return inputUri
    }

    static String normalizeQuotedCsvText(def value) {
        if (value == null) {
            return ""
        }

        // Convert to string, trim, and strip matching leading/trailing single or double quotes.
        return value.toString().trim().replaceAll(/^["']|["']$/, '')
    }

    static String deriveMobieViewName(def inputUri) {
        if (inputUri == null) {
            return ""
        }

        def text = inputUri.toString().trim()
        // Keep filename-like part even for directory-like paths/URIs.
        text = text.replaceAll(/\/+$/, '')
        if (!text) {
            return ""
        }

        def slashIndex = text.lastIndexOf('/')
        def viewName = slashIndex >= 0 ? text.substring(slashIndex + 1) : text

        // Drop URL query/fragment suffixes to keep names readable and stable.
        def queryIndex = viewName.indexOf('?')
        if (queryIndex >= 0) {
            viewName = viewName.substring(0, queryIndex)
        }
        def fragmentIndex = viewName.indexOf('#')
        if (fragmentIndex >= 0) {
            viewName = viewName.substring(0, fragmentIndex)
        }

        return viewName
    }

    static String shortStableHash(def value, int length = 8) {
        def text = value == null ? "" : value.toString()
        byte[] digest = MessageDigest.getInstance("SHA-1").digest(text.getBytes("UTF-8"))
        def hex = digest.collect { String.format("%02x", it) }.join('')
        return hex.substring(0, Math.min(length, hex.length()))
    }

    static createSegmentationMoBIECollectionTable(chInputRows, chLabels, outdir, outputFileName) {
        chInputRows
            .map { meta, image ->
                [meta.id, meta, image]
            }
            .join(
                chLabels.map { meta, label -> [meta.id, label] }
            )
            .flatMap { datasetID, meta, image, label ->
                def tableDir = Paths.get(outdir.toString()).toAbsolutePath().normalize()
                def inputUri = serializeMobieInputUri(image, tableDir)
                def baseViewName = deriveMobieViewName(inputUri)

                def labelRow = [
                    baseViewName,
                    "cellpose/${label.name}",
                    'OmeZarr',
                    '',
                    'labels',
                    'true',
                    ''
                ]

                def channelsText = normalizeQuotedCsvText(meta.segmentation_channels)
                def channelsForTable = channelsText ? channelsText.split(/\s*;\s*/).findAll { value -> value } : []

                def intensityRows
                if (channelsForTable) {
                    intensityRows = channelsForTable.collect { ch ->
                        [
                            baseViewName,
                            inputUri,
                            'OmeZarr',
                            ch.toString(),
                            'intensity',
                            'true',
                            'auto'
                        ]
                    }
                } else {
                    intensityRows = [
                        [
                            baseViewName,
                            inputUri,
                            'OmeZarr',
                            '',
                            'intensity',
                            'true',
                            'auto'
                        ]
                    ]
                }

                ([[sourceInputUri: inputUri, row: labelRow]] + intensityRows.collect { row -> [sourceInputUri: inputUri, row: row] })
            }
            .toList()
            .flatMap { annotatedRows ->
                def urisByBaseView = [:].withDefault { new LinkedHashSet<String>() }
                annotatedRows.each { annotated ->
                    def row = annotated.row
                    def sourceInputUri = annotated.sourceInputUri.toString()
                    urisByBaseView[row[0].toString()].add(sourceInputUri)
                }

                def resolvedViewByInputUri = [:]
                urisByBaseView.each { baseView, sourceUris ->
                    if (sourceUris.size() == 1) {
                        def onlyUri = sourceUris.iterator().next()
                        resolvedViewByInputUri[onlyUri] = baseView
                    } else {
                        sourceUris.each { sourceUri ->
                            resolvedViewByInputUri[sourceUri] = "${baseView}_${shortStableHash(sourceUri)}"
                        }
                    }
                }

                def seenIntensityRows = new LinkedHashSet<String>()
                annotatedRows.findAll { annotated ->
                    def row = annotated.row
                    def sourceInputUri = annotated.sourceInputUri.toString()
                    row[0] = resolvedViewByInputUri[sourceInputUri]

                    def type = row[4]
                    if (type != 'intensity') {
                        return true
                    }

                    // Same input URI rows should share a view, so keep each intensity source once per view+channel.
                    def dedupKey = row.join('\t')
                    return seenIntensityRows.add(dedupKey)
                }.collect { annotated -> annotated.row }
            }
            .map { view, uri, format, channel, type, exclusive, contrastLimits ->
                "${view},${uri},${format},${channel},${type},${exclusive},${contrastLimits}"
            }
            .collectFile(
                name: outputFileName,
                storeDir: outdir,
                seed: 'view,uri,format,channel,type,exclusive,contrast_limits',
                newLine: true
            )
    }

    static String normalizeOptionalText(def value) {
        if (value == null) {
            return null
        }

        def text = value.toString().trim()
        text = text.replaceAll(/^"+|"+$/, '')
        text = text.replaceAll(/^'+|'+$/, '')
        text = text.trim()

        return text ? text : null
    }

    static String optionalCliArg(String flag, def value) {
        def normalized = normalizeOptionalText(value)
        return normalized != null ? "${flag} ${normalized}" : ""
    }

    static String requireColumnValue(def row, String columnName) {
        def value = row[columnName]
        if (value == null || !value.toString().trim()) {
            throw new IllegalArgumentException("Samplesheet is missing required '${columnName}' column value.")
        }
        def text = value.toString().trim()
        // Tolerate accidental quote artifacts in CSV fields.
        text = text.replaceAll(/^"+/, '').replaceAll(/"+$/, '')
        text = text.replaceAll(/^'+/, '').replaceAll(/'+$/, '')
        return text
    }

    static String resolvePath(String uri, String samplesheetPath) {
        
        if (uri == null) {
            return null
        }

        // absolute uri
        def hasUriScheme = (uri =~ /(?i)https?:/)
        if(hasUriScheme) {
            return uri
        }

        // absolute file path
        Path inputPath = Paths.get(uri)
        if(inputPath.isAbsolute()) {
            if (!Files.exists(inputPath)) {
                throw new IllegalArgumentException("No such file or directory: ${inputPath}")
            }    
            return inputPath
        }

        // relative file path -- resolve to absolute so downstream `path()` process inputs
        // (which need a properly resolvable filesystem path, not one relative to whatever
        // the task's working directory happens to be) can stage it correctly.
        Path resolvedPath = Paths.get(samplesheetPath).toAbsolutePath().parent.resolve(inputPath).normalize()
        if (!Files.exists(resolvedPath)) {
            throw new IllegalArgumentException("No such file or directory: ${resolvedPath}")
        }
        return resolvedPath
    }

    static List<Integer> parseIndexList(def value, String columnName) {
        if (!value) {
            return null
        }

        def text = value.toString().trim()
        text = text.replaceAll(/^"+|"+$/, '')
        text = text.replaceAll(/^'+|'+$/, '')
        if (!text) {
            return null
        }

        def tokens = text.split(/[;,\s]+/).findAll { it }
        if (!tokens) {
            return null
        }

        try {
            return tokens.collect { token -> Integer.parseInt(token) }
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException(
                "Invalid '${columnName}' value '${value}'. Expected integer indices separated by comma, semicolon, or whitespace.",
                e
            )
        }
    }

    static List<Integer> parseSemicolonIndexList(def value, String columnName) {
        if (!value) {
            return null
        }

        def text = value.toString().trim()
        text = text.replaceAll(/^"+|"+$/, '')
        text = text.replaceAll(/^'+|'+$/, '')
        if (!text) {
            return null
        }

        def tokens = text.split(/\s*;\s*/).findAll { it }
        if (!tokens) {
            return null
        }

        try {
            return tokens.collect { token -> Integer.parseInt(token) }
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException(
                "Invalid '${columnName}' value '${value}'. Expected integer indices separated by semicolon (for example: 0;1;2).",
                e
            )
        }
    }

    static List<Double> parseFloatList(def value, String columnName, Integer expectedCount = null) {
        if (!value) {
            return null
        }

        def text = value.toString().trim()
        text = text.replaceAll(/^"+|"+$/, '')
        text = text.replaceAll(/^'+|'+$/, '')
        if (!text) {
            return null
        }

        def tokens = text.split(/[;,\s]+/).findAll { it }
        if (!tokens) {
            return null
        }

        if (expectedCount != null && tokens.size() != expectedCount) {
            throw new IllegalArgumentException(
                "Invalid '${columnName}' value '${value}'. Expected exactly ${expectedCount} numeric values."
            )
        }

        try {
            return tokens.collect { token -> Double.parseDouble(token) }
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException(
                "Invalid '${columnName}' value '${value}'. Expected numeric values separated by comma, semicolon, or whitespace.",
                e
            )
        }
    }

    static Integer parseOptionalInt(def value, String columnName) {
        if (!value) {
            return null
        }

        def text = value.toString().trim()
        text = text.replaceAll(/^"+|"+$/, '')
        text = text.replaceAll(/^'+|'+$/, '')
        if (!text) {
            return null
        }

        try {
            return Integer.parseInt(text)
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException(
                "Invalid '${columnName}' value '${value}'. Expected an integer.",
                e
            )
        }
    }

    static Double parseOptionalFloat(def value, String columnName) {
        if (!value) {
            return null
        }

        def text = value.toString().trim()
        text = text.replaceAll(/^"+|"+$/, '')
        text = text.replaceAll(/^'+|'+$/, '')
        if (!text) {
            return null
        }

        try {
            return Double.parseDouble(text)
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException(
                "Invalid '${columnName}' value '${value}'. Expected a numeric value.",
                e
            )
        }
    }

    static List<Double> parseOptionalFloatRange(def value, String columnName) {
        if (!value) {
            return null
        }

        def text = value.toString().trim()
        text = text.replaceAll(/^"+|"+$/, '')
        text = text.replaceAll(/^'+|'+$/, '')
        if (!text) {
            return null
        }

        def matcher = (text =~ /^\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*:\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*$/)
        if (!matcher.matches()) {
            throw new IllegalArgumentException(
                "Invalid '${columnName}' value '${value}'. Expected 'start:end' with numeric bounds (for example: 12.5:48.0)."
            )
        }

        def start = Double.parseDouble(matcher[0][1])
        def end = Double.parseDouble(matcher[0][2])
        if (end <= start) {
            throw new IllegalArgumentException(
                "Invalid '${columnName}' value '${value}'. End value must be greater than start value."
            )
        }

        return [start, end]
    }

    static List<Integer> parseOptionalRange(def value, String columnName) {
        if (!value) {
            return null
        }

        def text = value.toString().trim()
        text = text.replaceAll(/^"+|"+$/, '')
        text = text.replaceAll(/^'+|'+$/, '')
        if (!text) {
            return null
        }

        def matcher = (text =~ /^\s*(\d+)\s*:\s*(\d+)\s*$/)
        if (!matcher.matches()) {
            throw new IllegalArgumentException(
                "Invalid '${columnName}' value '${value}'. Expected 'start:end' with integer bounds (for example: 128:512)."
            )
        }

        def start = Integer.parseInt(matcher[0][1])
        def end = Integer.parseInt(matcher[0][2])
        if (start < 0) {
            throw new IllegalArgumentException(
                "Invalid '${columnName}' value '${value}'. Start index must be >= 0."
            )
        }
        if (end <= start) {
            throw new IllegalArgumentException(
                "Invalid '${columnName}' value '${value}'. End index must be greater than start index."
            )
        }

        return [start, end]
    }
}


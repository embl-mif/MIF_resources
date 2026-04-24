import argparse
import dask.array as da
import pandas as pd
from skimage.measure import regionprops_table
import pyarrow as pa
import pyarrow.parquet as pq

def run_measurement(input_zarr, output_parquet, intensity_path, label_path):
    dask.config.set(scheduler='threads')
    intensity_stack = da.from_zarr(input_zarr, component=intensity_path)
    label_stack = da.from_zarr(input_zarr, component=label_path)
    
    n_tps = label_stack.shape[0]
    features = ["label", "area", "centroid", "intensity_mean", "intensity_max"]
    
    tables = []
    
    for t in range(n_tps):
        print(f"Processing timepoint {t+1}/{n_tps}...")
        labels_t = label_stack[t].compute()
        intensity_t = intensity_stack[t].compute()

        # Transform
        intensity_t_moved = da.moveaxis(intensity_t, 0, -1).compute()

        # Measure
        table = regionprops_table(labels_t, intensity_t_moved, properties=features)
        
        # Bookkeeping
        df = pd.DataFrame(table)
        df["t"] = t 
        tables.append(df)
        
    final_df = pd.concat(tables, ignore_index=True)
    
    # Convert DataFrame to PyArrow Table to ensure explicit types
    table = pa.Table.from_pandas(final_df)
    pq.write_table(table, output_parquet, compression='snappy')
    
    print(f"Successfully saved to {output_parquet}")
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract features from OME-Zarr for Nextflow.")
    parser.add_argument("-z", "--input_zarr", type=str, required=True, help="Input Zarr path")
    parser.add_argument("-o", "--output_path", type=str, required=True, help="Output Parquet path")
    parser.add_argument("--intensity", type=str, default="0", help="Component for intensity (e.g., '0')")
    parser.add_argument("--labels", type=str, default="labels/nuclei", help="Component for labels (e.g., 'labels/nuclei')")

    args = parser.parse_args()
    run_measurement(args.input_zarr, args.output_path, args.intensity, args.labels)

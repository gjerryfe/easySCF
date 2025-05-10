#!/usr/bin/env python3
import argparse
import os
import sys
import anndata
import h5py

def convert_h5_to_h5ad(input_path, output_path):
    """Convert H5 file to h5ad format."""
    try:
        with h5py.File(input_path, 'r') as h5_file:
            # Create AnnData object with all components
            adata = anndata.AnnData(
                X=h5_file['X'][:],
                obs=dict(h5_file['obs']),
                var=dict(h5_file['var']),
                uns=dict(h5_file['uns']),
                obsm={k: v[:] for k, v in h5_file['obsm'].items()},
                varm={k: v[:] for k, v in h5_file['varm'].items()}
            )
        adata.write_h5ad(output_path)
        print(f"Successfully converted {input_path} to {output_path}")
    except Exception as e:
        print(f"Error converting H5 to h5ad: {str(e)}", file=sys.stderr)
        sys.exit(1)

def convert_h5ad_to_h5(input_path, output_path):
    """Convert h5ad file to H5 format."""
    try:
        adata = anndata.read_h5ad(input_path)
        with h5py.File(output_path, 'w') as h5_file:
            # Store all AnnData components
            h5_file.create_dataset('X', data=adata.X)
            
            # Store obs, var, uns
            for key, value in adata.obs.items():
                h5_file.create_dataset(f'obs/{key}', data=value.values)
            for key, value in adata.var.items():
                h5_file.create_dataset(f'var/{key}', data=value.values)
            for key, value in adata.uns.items():
                h5_file.create_dataset(f'uns/{key}', data=str(value))
                
            # Store obsm and varm
            for key, value in adata.obsm.items():
                h5_file.create_dataset(f'obsm/{key}', data=value)
            for key, value in adata.varm.items():
                h5_file.create_dataset(f'varm/{key}', data=value)
                
        print(f"Successfully converted {input_path} to {output_path}")
    except Exception as e:
        print(f"Error converting h5ad to H5: {str(e)}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="Convert between H5 and h5ad (AnnData) formats"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input file (.h5 or .h5ad)"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output file path (extension determines format)"
    )
    
    args = parser.parse_args()
    
    # Validate input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} does not exist!", file=sys.stderr)
        sys.exit(1)
    
    # Determine conversion direction based on file extensions
    input_ext = os.path.splitext(args.input)[1].lower()
    output_ext = os.path.splitext(args.output)[1].lower()
    
    if input_ext == '.h5' and output_ext == '.h5ad':
        convert_h5_to_h5ad(args.input, args.output)
    elif input_ext == '.h5ad' and output_ext == '.h5':
        convert_h5ad_to_h5(args.input, args.output)
    else:
        print("Error: Input/output extensions must be .h5/.h5ad pair", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()

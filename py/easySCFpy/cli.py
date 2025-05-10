#!/usr/bin/env python3
import argparse
import os
import sys
from pathlib import Path
from easySCFpy import loadH5, saveH5
import anndata

def convert_h5_to_h5ad(input_path, output_path):
    """Convert H5 file to h5ad format."""
    try:
        adata = loadH5(input_path)
        adata.write_h5ad(output_path)
        print(f"Successfully converted {input_path} to {output_path}")
    except Exception as e:
        print(f"Error converting H5 to h5ad: {str(e)}", file=sys.stderr)
        sys.exit(1)

def convert_h5ad_to_h5(input_path, output_path):
    """Convert h5ad file to H5 format."""
    try:
        adata = anndata.read_h5ad(input_path)
        saveH5(adata, output_path)
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

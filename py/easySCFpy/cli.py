#!/usr/bin/env python3
import argparse
import os
import sys
import subprocess
from pathlib import Path
from easySCFpy import loadH5, saveH5
import anndata

def detect_file_type(file_path):
    """Detect if file is R or Python format"""
    ext = os.path.splitext(file_path)[1].lower()
    if ext == '.rds':
        return 'R'
    elif ext in ('.h5', '.h5ad'):
        return 'Python'
    else:
        raise ValueError(f"Unsupported file type: {ext}")

def call_r_script(input_path, output_path, operation):
    """Call R script for conversion"""
    # r_script = os.path.join(os.path.dirname(__file__), '../r/inst/exec/easy-scf-tool-r')
    cmd = ['easy-scf-tool-r', '-i', input_path, '-o', output_path]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"R script failed: {result.stderr}", file=sys.stderr)
        sys.exit(1)
    print(result.stdout)

def convert_h5_to_h5ad(input_path, output_path):
    """Convert H5 file to h5ad format using Python implementation"""
    try:
        adata = loadH5(input_path)
        adata.write_h5ad(output_path)
        print(f"Successfully converted {input_path} to {output_path}")
    except Exception as e:
        print(f"Error converting H5 to h5ad: {str(e)}", file=sys.stderr)
        sys.exit(1)

def convert_h5ad_to_h5(input_path, output_path):
    """Convert h5ad file to H5 format using Python implementation"""
    try:
        adata = anndata.read_h5ad(input_path)
        saveH5(adata, output_path)
        print(f"Successfully converted {input_path} to {output_path}")
    except Exception as e:
        print(f"Error converting h5ad to H5: {str(e)}", file=sys.stderr)
        sys.exit(1)

def convert_rds_to_h5(input_path, output_path):
    """Convert RDS file to H5 format using R implementation"""
    call_r_script(input_path, output_path, 'rds_to_h5')

def convert_h5_to_rds(input_path, output_path):
    """Convert H5 file to RDS format using R implementation"""
    call_r_script(input_path, output_path, 'h5_to_rds')

def main():
    parser = argparse.ArgumentParser(
        description="Unified converter between R and Python single-cell formats"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input file (.rds, .h5 or .h5ad)"
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
    
    # Python conversions
    if input_ext == '.h5' and output_ext == '.h5ad':
        convert_h5_to_h5ad(args.input, args.output)
    elif input_ext == '.h5ad' and output_ext == '.h5':
        convert_h5ad_to_h5(args.input, args.output)
    # R conversions
    elif input_ext == '.rds' and output_ext == '.h5':
        convert_rds_to_h5(args.input, args.output)
    elif input_ext == '.h5' and output_ext == '.rds':
        convert_h5_to_rds(args.input, args.output)
    # R+Python combined conversion (RDS -> H5 -> h5ad)
    elif input_ext == '.rds' and output_ext == '.h5ad':
        import tempfile
        with tempfile.NamedTemporaryFile(suffix='.h5') as tmp:
            convert_rds_to_h5(args.input, tmp.name)
            convert_h5_to_h5ad(tmp.name, args.output)
    else:
        print("Error: Unsupported conversion pair", file=sys.stderr)
        print("Supported conversions:", file=sys.stderr)
        print("  Python: .h5 <-> .h5ad", file=sys.stderr)
        print("  R: .rds <-> .h5", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()

#!/bin/bash

python -m iedb_immrep \
    --export_file /home/lwoods/workspace/IEDB_IMMREP/dat/api_export_results_HUMAN_I.parquet \
    --output_dir dat/HUMAN_I \
    --species human

python -m iedb_immrep \
    --export_file /home/lwoods/workspace/IEDB_IMMREP/dat/api_export_results_HUMAN_II.parquet \
    --output_dir dat/HUMAN_II \
    --species human

# python -m iedb_immrep \
#     --export_file /home/lwoods/workspace/IEDB_IMMREP/dat/api_export_results_MOUSE_I.parquet \
#     --output_dir dat/MOUSE_I \
#     --species mouse

# python -m iedb_immrep \
#     --export_file /home/lwoods/workspace/IEDB_IMMREP/dat/api_export_results_MOUSE_II.parquet \
#     --output_dir dat/MOUSE_II \
#     --species mouse

# python iedb_immrep/process_vdjb.py
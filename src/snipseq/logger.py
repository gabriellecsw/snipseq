"""
@author      Gabrielle C.
@create date 2026-03-17
@desc        Logging functions for Snipseq
"""

import logging
import typer

def log(msg: str):
    typer.echo(msg)
    logging.info(msg)

def log_error(msg: str):
    typer.echo(msg, err=True)
    logging.error(msg)

def setup_logger(log_file):
    logging.basicConfig(filename=log_file, 
        level=logging.INFO, 
        format="%(asctime)s - %(message)s", 
        datefmt="%Y-%m-%d %H:%M:%S")
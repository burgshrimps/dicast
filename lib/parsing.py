import argparse
import sys

def parse_arguments(arguments = sys.argv[1:]):
    parser = argparse.ArgumentParser(description='dicast')

    subparsers = parser.add_subparsers(dest='command', help='Subcommands')

    parser_prepare = subparsers.add_parser('prepare', help='Variant Preparation and Feature Extraction.')
    parser_prepare.add_argument('sample', help='Sample name')
    parser_prepare.add_argument('ref', help='Reference genome name', default='hg38')
    parser_prepare.add_argument('params', help='Preparation parameter file')
    parser_prepare.add_argument('workdir', help='Working and output directory')

    parser_train = subparsers.add_parser('train', help='Train a new Model.')
    parser_train.add_argument('params', help='Training parameter file')


    parser_predict = subparsers.add_parser('predict', help='Make Predictions on a Set of Variants..')

    return parser.parse_args(arguments)
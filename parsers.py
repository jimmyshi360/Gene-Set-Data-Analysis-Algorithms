
class Parsers:

    def __init__(self,parsers):
        from argparse import ArgumentParser

        usage = "usage: %prog [options]"
        parser = ArgumentParser(usage, version="%prog dev-unreleased")

        if parsers.count("-g") == 1:
            parser.add_argument(
                "-g",
                "--gene-list file",
                dest="gene_list",
                help="gene list file for comparision",
                metavar=".gmt FILE",
                required=True
            )
        if parsers.count("-a") == 1:
            parser.add_argument(
                "-a",
                "--annotations-file",
                dest="annotation_list",
                help="annotation file",
                metavar=".gmt FILE",
                required=True
            )
        if parsers.count("-b") == 1:
            parser.add_argument(
                "-b",
                "--background-gene file",
                dest="background_list",
                help="background gene list file for comparision (OPTIONAL)",
                metavar=".txt FILE"
            )
        if parsers.count("-c") == 1:
            parser.add_argument(
                "-c",
                "--clusters-file",
                dest="cluster_list",
                help="cluster file",
                metavar=".mat FILE",
                required=True
            )
        if parsers.count("-l") == 1:
            parser.add_argument(
                "-l",
                "--cluster number",
                dest="cluster_number",
                help="the cluster to run through (DEFAULT 0)",
                metavar="INTEGER",
                type=int,
                default=0)
        if parsers.count("-o") == 1:
            parser.add_argument(
                "-o",
                "--output-file",
                dest="output",
                help="file to output",
                metavar=".txt FILE",
                required=True
            )
        if parsers.count("-r") == 1:
            parser.add_argument(
                "-r",
                "--the alpha level for the test",
                dest="rate",
                help="a decimal for the alpha rate (DEFAULT 0.05)",
                metavar="FLOAT",
                type=float,
                default=0.05)
        self._args = parser.parse_args()

    @property
    def args(self):
        return self._args
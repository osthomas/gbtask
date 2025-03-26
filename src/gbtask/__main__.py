import logging

from gbtask.cli import parse_args


def main():
    args = parse_args()
    levels = sorted(set(logging.getLevelNamesMapping().values()))
    loglevel = levels[-min(len(levels), args.v + 1)]
    logging.basicConfig(
        level=loglevel, format="%(levelname)s\t%(name)-20s:\t%(message)s"
    )
    log = logging.getLogger(__name__)
    args.func(args)


if __name__ == "__main__":
    main()

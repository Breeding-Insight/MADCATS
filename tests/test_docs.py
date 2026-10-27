"""Keep the command reference pages in step with the command-line parser.

Every option a subcommand accepts must appear in a table on its page in
docs/src/commands/, and every option those tables list must exist. Where a
table gives a choice set or a plain default, it must match the parser too.
"""
import argparse
import re
import unittest
from pathlib import Path

import madcats

COMMANDS_DIR = Path(__file__).resolve().parents[1] / "docs" / "src" / "commands"
OPTION = re.compile(r"`(--[a-z0-9-]+)")
CHOICES = re.compile(r"\{([^}]*)\}")


def subcommand_parsers():
    parser = madcats.build_parser()
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            return action.choices
    raise AssertionError("madcats has no subcommands")


def documented_rows(page: Path):
    """Yield (options, first cell, default cell or None) for each table row naming options."""
    for line in page.read_text().splitlines():
        if not line.startswith("| `--"):
            continue
        cells = [cell.strip() for cell in re.split(r"(?<!\\)\|", line)[1:-1]]
        options = OPTION.findall(cells[0])
        default = cells[1] if len(cells) == 3 else None
        yield options, cells[0], default


def documented_default_matches(action: argparse.Action, cell: str):
    """True or False when the cell can be compared with the parser default, None when it cannot."""
    text = cell.strip("`")
    if text == "none":
        return action.default in (None, "")
    if isinstance(action.default, bool):
        return text == ("on" if action.default else "off")
    if isinstance(action.default, (int, float)):
        try:
            return float(text) == float(action.default)
        except ValueError:
            return None  # descriptive defaults such as "up to 4"
    if isinstance(action.default, str) and " " not in text:
        return text == action.default
    return None


@unittest.skipUnless(COMMANDS_DIR.is_dir(), "documentation sources not present")
class CommandReferenceTests(unittest.TestCase):
    def test_option_tables_match_parser(self):
        for name, subparser in subcommand_parsers().items():
            with self.subTest(command=name):
                page = COMMANDS_DIR / f"{name}.md"
                self.assertTrue(page.is_file(), f"no reference page for {name}")
                actions = {
                    option: action
                    for action in subparser._actions
                    for option in action.option_strings
                    if option.startswith("--") and option != "--help"
                }
                rows = list(documented_rows(page))
                documented = {option for options, _, _ in rows for option in options}
                self.assertEqual(sorted(set(actions) - documented), [], f"options missing from {page.name}")
                self.assertEqual(sorted(documented - set(actions)), [], f"options in {page.name} that {name} lacks")

                for options, first_cell, default in rows:
                    if len(options) != 1:
                        continue
                    action = actions[options[0]]
                    listed = CHOICES.search(first_cell)
                    if listed and action.choices is not None:
                        self.assertEqual(
                            listed.group(1).split(","),
                            [str(choice) for choice in action.choices],
                            f"{options[0]} choices in {page.name}",
                        )
                    if default is not None:
                        self.assertIsNot(
                            documented_default_matches(action, default),
                            False,
                            f"{options[0]} default in {page.name}: docs say {default!r}, parser has {action.default!r}",
                        )


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3

# Python standard library imports
import argparse
import glob
import re
import sys


# Rules that are applied before stripping away comments
_pre_rules = [
    {
        "name": "comment_space",
        "description": "There shall be a space after the comment marker",
        "pattern": r"^( *)(!(?!\$))([^ \n])",
        "replace": r"\1\2 \3",
    },
    {
        "name": "line_continuation",
        "description": "Lines shall not end with \\",
        "pattern": r"\\ *$",
        "replace": None,
    },
    {
        "name": "trailing_whitespace",
        "description": "Lines shall not have trailing whitespace",
        "pattern": r"[ \t]+$",
        "replace": "",
    },
]


# Rules that apply only to code that are not comments
_rules = [
    {
        "name": "space_after_comma",
        "description": "Comma should be followed by a whitespace",
        "pattern": r"(,)([^ ])",
        "replace": r", \2",
    },
    {
        "name": "space_before_comma",
        "description": "No whitespace before comma",
        "pattern": r" +,",
        "replace": r",",
    },
    {
        "name": "end_subroutine_1",
        "description": "Use END SUBROUTINE instead of ENDSUBROUTINE",
        "pattern": r"^( *)(ENDSUBROUTINE)",
        "replace": r"\1END SUBROUTINE",
    },
    {
        "name": "end_subroutine_2",
        "description": "END SUBROUTINE shall be followed by subroutine name",
        "pattern": r"^ *END *SUBROUTINE *$",
        "replace": None,
    },
    {
        "name": "end_function_1",
        "description": "Use END FUNCTION instead of ENDFUNCTION",
        "pattern": r"^( *)(ENDFUNCTION)",
        "replace": r"\1END FUNCTION",
    },
    {
        "name": "end_function_2",
        "description": "END FUNCTION shall be followed by function name",
        "pattern": r"^ *END *FUNCTION *$",
        "replace": None,
    },
    {
        "name": "end_do",
        "description": "Use END DO instead of ENDDO",
        "pattern": r"^( *)(ENDDO)",
        "replace": r"\1END DO",
    },
    {
        "name": "end_if",
        "description": "Use END IF instead of ENDIF",
        "pattern": r"^( *)(ENDIF)",
        "replace": r"\1END IF",
    },
    {
        "name": "end",
        "description": "Use END <construct> instead of END",
        "pattern": r"^ *END *$",
        "replace": None,
    },
    {
        "name": "open_parenthesis",
        "description": "Opening parenthesis shall not be followed by a space",
        "pattern": r"\( +(?!&)",
        "replace": "(",
    },
    {
        "name": "close_parenthesis",
        "description": "Closing parenthesis shall not follow a space",
        "pattern": r" +\)",
        "replace": ")",
    },
    {
        "name": "tab",
        "description": "Tab is not allowed in Fortran code",
        "pattern": r"\t",
        "replace": None,
    },
    {
        "name": "eq",
        "description": "Old-style Fortran operator shall not be used",
        "pattern": r"\.eq\.",
        "replace": "==",
    },
    {
        "name": "neq",
        "description": "Old-style Fortran operator shall not be used",
        "pattern": r"\.neq\.",
        "replace": "/=",
    },
    {
        "name": "leq",
        "description": "Old-style Fortran operator shall not be used",
        "pattern": r"\.leq\.",
        "replace": "<=",
    },
    {
        "name": "geq",
        "description": "Old-style Fortran operator shall not be used",
        "pattern": r"\.geq\.",
        "replace": ">=",
    },
    {
        "name": "lt",
        "description": "Old-style Fortran operator shall not be used",
        "pattern": r"\.lt\.",
        "replace": "<",
    },
    {
        "name": "gt",
        "description": "Old-style Fortran operator shall not be used",
        "pattern": r"\.gt\.",
        "replace": ">",
    },
]


def comment_pos(line):
    """ Return the position of the comment marker, taking strings into
    account. In Fortran, we cannot escape quotation marks or apostrophes.
    Some compilers implement this as extensions, but it is not standard and
    should never be used."""

    # We need to keep track of whether we are inside a string or not.
    stack = []

    # Walk string. When a " or ' is found, push to stack. When the same
    # character is found, pop from stack. If stack is empty, we are not in a
    # string.
    pos = len(line)
    for i, c in enumerate(line):
        if c in ["'", '"']:
            if stack and stack[-1] == c:
                stack.pop()
            else:
                stack.append(c)
        elif c == "!" and not stack:
            pos = i

    # Wind back pos when there are whitespaces before the comment marker
    while pos > 0 and line[pos-1] in [" ", "\t"]:
        pos -= 1

    return pos


def lint_line(fname, lineno, line):
    pass_one = True

    # First pass is the rules we apply before stripping comments
    for rule in _pre_rules:
        pattern = rule["pattern"]
        if re.search(pattern, line):
            pass_one = False
            print(f"::error file={fname},line={lineno}::{rule['description']}")

    # Identify position of comment marker ! taking strings into account
    # Example of a problematic line:
    #     WRITE(*, '("Attribute ", A, " does not exist!")') key
    pos = comment_pos(line)
    line_wo_comments = line[:pos]

    # Second pass is the rules we apply to the code after stripping comments
    for rule in _rules:
        pattern = rule["pattern"]
        if re.search(pattern, line_wo_comments):
            pass_one = False
            print(f"::error file={fname},line={lineno}::{rule['description']}")

    # We allow lines to be 80 characters long - counting only code. Comments
    # are not counted. Therefore we need a special function to count the line
    # length.
    maxlen = 80
    thislen = len(line_wo_comments.rstrip("\n"))
    if thislen > maxlen:
        print(f"::error file={fname},line={lineno}"
              f"::Line too long {thislen} > {maxlen}")
        pass_one = False

    return pass_one


def get_patch():
    import whatthepatch

    """ Read patch from stdin """
    diff = sys.stdin.read()
    return whatthepatch.parse_patch(diff)


def lint_patch(patch):
    pass_all = True
    for diff in patch:
        fname = diff.header.new_path

        # When the file is not a Fortran file, skip it
        if not fname.endswith(".F90"):
            continue

        for change in diff.changes:
            # No need to check for removed lines
            if change.new is None:
                continue

            # Check that line
            pass_one = lint_line(fname, change.new, change.line)
            if not pass_one:
                pass_all = False

    return pass_all


def lint_file(fname):
    if not fname.endswith(".F90"):
        raise ValueError(f"{fname} is not a Fortran file")

    pass_all = True
    with open(fname, "r") as f:
        for no, line in enumerate(f):
            pass_one = lint_line(fname, no+1, line)
            if not pass_one:
                pass_all = False

    return pass_all


def fix_line(line):
    """ Fix a Fortran line - match regex and replace """
    # First apply rules that are applied before stripping comments
    for rule in _pre_rules:
        pattern = rule["pattern"]
        replace = rule["replace"]
        if re.search(pattern, line):
            if replace is not None:
                line = re.sub(pattern, replace, line)

    # Then apply rules that apply to the code after stripping comments
    for rule in _rules:
        # The position of the comment marker can change as rules are applied
        pos = comment_pos(line)
        pattern = rule["pattern"]
        replace = rule["replace"]
        if re.search(pattern, line[:pos]):
            if replace is not None:
                line = re.sub(pattern, replace, line[:pos]) + line[pos:]

    return line


def fix_file(fname):
    """ Fix a Fortran file - match regex and replace """
    if not fname.endswith(".F90"):
        raise ValueError(f"{fname} is not a Fortran file")

    with open(fname, "r") as f:
        lines = f.readlines()

    with open(fname, "w") as f:
        for line in lines:
            f.write(fix_line(line))

    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stdin", action="store_true")
    parser.add_argument("--file", type=str, default=None)
    parser.add_argument("--fix", type=str, default=None)
    parser.add_argument("--fixall", action="store_true")
    args = parser.parse_args()

    passed_all = True
    if args.stdin:
        patch = get_patch()
        passed_all = lint_patch(patch)
    elif args.file:
        passed_all = lint_file(args.file)
    elif args.fix:
        fix_file(args.fix)
    elif args.fixall:
        for filename in glob.iglob("**/*.F90", recursive=True):
            fix_file(filename)
    else:
        for filename in glob.iglob("**/*.F90", recursive=True):
            if not lint_file(filename):
                passed_all = False

    sys.exit(0 if passed_all else 1)


if __name__ == "__main__":
    main()

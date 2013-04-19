#!/usr/bin/env python
"""
    Subversion pre-commit hook for sand.
    Checks that:
    * source files does not contain tabs.
    * there are no files ending with whitespace in source files.

    Usage: %prog REPOS TXN

    If the -r or --revision option is given TXN refers to a revision.
    This is used for testing.
"""

import argparse
import os
import pipes # For pipes.quote. In 3.2 we can use shlex.quote instead.
import re
import shlex
import subprocess

# file extensions for code files
code_extensions = ".cpp .hpp .c .h *.C *.H .cu .py .tcc .cc"
modified = re.compile('^(?:M|A)(\s+)(?P<name>.*)')

def command_output(cmd):
    " Capture a command's standard output. "
    return subprocess.Popen(
        shlex.split(cmd), stdout=subprocess.PIPE).communicate()[0]


def files_changed(all_files):
    """ List the files added or updated by this transaction.

    If all_files is true, all files in repository are returned.
    """
    files = []
    if all_files:
        for root, dirs, file_names in os.walk('.'):
            for file_name in file_names:
                files.append(os.path.join(root, file_name))
    else:
        p = subprocess.Popen(['git', 'status', '--porcelain'], stdout=subprocess.PIPE)
        out, err = p.communicate()
        for line in out.splitlines():
            match = modified.match(line)
            if match:
                files.append(match.group('name'))
    return files


def file_contents(filename):
    " Return a file's contents for this transaction. "
    return command_output("cat %s" % (pipes.quote(filename)))


def file_diff(filename):
    "Return the changed lines for the given file. "
    diff = command_output("git diff --cached -U0 %s" % (pipes.quote(filename)))
    new_lines = ''
    for line in diff.split('\n'):
        if line.startswith('+'):
            new_lines += line + '\n'
    return new_lines


def contains_tabs(file_content):
    " Return True if this version of the file contains tabs. "
    return "\t" in file_content


def contains_trailing_whitespace(file_content):
    " Return True if this version of the file contains lines ending with a whitespace. "
    for line in file_content.splitlines():
        if line.endswith(' '): # One space is enough, tab is checked elsewhere
            return True
    return False


def check_code_files(all_files):
    " Check C++ files in this transaction are tab-free. "
    def is_code_file(fname):
        import os
        return os.path.splitext(fname)[1] in code_extensions.split()
    files_with_tabs = []
    files_with_trailing_whitespace = []
    for ff in files_changed(all_files):
        if is_code_file(ff):
            content_to_check = ''
            if all_files:
                content_to_check = file_contents(ff)
            else:
                content_to_check = file_diff(ff)

            if (contains_tabs(content_to_check)):
                files_with_tabs.append(ff)
            if (contains_trailing_whitespace(content_to_check)):
                files_with_trailing_whitespace.append(ff)
    if len(files_with_tabs) > 0:
        sys.stderr.write("The following files contain tabs:\n%s\n\n"
                          % "\n".join(files_with_tabs))
    if len(files_with_trailing_whitespace) > 0:
        sys.stderr.write("The following files contains lines with trailing whitespace:\n%s\n\n"
                          % "\n".join(files_with_trailing_whitespace))
        sys.stderr.write('Hint: Reg Ex "[\\t ]+$" - to find trailing white spaces\n\n')
    return len(files_with_tabs) + len(files_with_trailing_whitespace)


def main():
    parser = argparse.ArgumentParser(description='Pre-commit hook for git.')
    parser.add_argument('-a', '--all-files',
                        help="Test all files in repository, not just the changed ones.",
                        action="store_true", default=False)
    errors = 0
    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        errors += 1
    else:
        # No need to stash, since we only use 'git diff --cached'.
        errors += check_code_files(args.all_files)

    return errors


if __name__ == "__main__":
  import sys
  sys.exit(main())


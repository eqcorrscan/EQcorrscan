"""
Helper functions for access to archives, e.g. Party, Tribe

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""

import contextlib
import os
import shutil
import tempfile
import tarfile
import logging


Logger = logging.getLogger(__name__)


@contextlib.contextmanager
def temporary_directory():
    """ make a temporary directory, yeild its name, cleanup on exit """
    dir_name = tempfile.mkdtemp()
    yield dir_name
    if os.path.exists(dir_name):
        shutil.rmtree(dir_name)


def _par_read(dirname, compressed=True):
    """
    Internal write function to read a formatted parameter file.

    :type dirname: str
    :param dirname: Directory to read the parameter file from.
    :type compressed: bool
    :param compressed: Whether the directory is compressed or not.
    """
    from eqcorrscan.core.match_filter.matched_filter import MatchFilterError
    from eqcorrscan.core.match_filter.template import Template

    templates = []
    if compressed:
        arc = tarfile.open(dirname, "r:*")
        members = arc.getmembers()
        _parfile = [member for member in members
                    if member.name.split(os.sep)[-1] ==
                    'template_parameters.csv']
        if len(_parfile) == 0:
            arc.close()
            raise MatchFilterError(
                'No template parameter file in archive')
        parfile = arc.extractfile(_parfile[0])
    else:
        parfile = open(dirname + '/' + 'template_parameters.csv', 'r')
    for line in parfile:
        t_in = Template()
        for key_pair in line.rstrip().split(','):
            if key_pair.split(':')[0].strip() == 'name':
                t_in.__dict__[key_pair.split(':')[0].strip()] = \
                    key_pair.split(':')[-1].strip()
            elif key_pair.split(':')[0].strip() == 'filt_order':
                try:
                    t_in.__dict__[key_pair.split(':')[0].strip()] = \
                        int(key_pair.split(':')[-1])
                except ValueError:
                    pass
            else:
                try:
                    t_in.__dict__[key_pair.split(':')[0].strip()] = \
                        float(key_pair.split(':')[-1])
                except ValueError:
                    pass
        templates.append(t_in)
    parfile.close()
    if compressed:
        arc.close()
    return templates


def _resolved(x):
    return os.path.realpath(os.path.abspath(x))


def _badpath(path, base):
    """
    joinpath will ignore base if path is absolute.
    """
    return not _resolved(os.path.join(base, path)).startswith(base)


def _badlink(info, base):
    """
    Links are interpreted relative to the directory containing the link
    """
    tip = _resolved(os.path.join(base, os.path.dirname(info.name)))
    return _badpath(info.linkname, base=tip)


def _safemembers(members):
    """
    Check members of a tar archive for safety.

    Ensure that they do not contain paths or links outside of where we
    need them - this would only happen if the archive wasn't made by
    eqcorrscan.

    :type members: :class:`tarfile.TarFile`
    :param members: an open tarfile.
    """
    base = _resolved(".")

    for finfo in members:
        if _badpath(finfo.name, base):
            print(finfo.name, "is blocked (illegal path)")
        elif finfo.issym() and _badlink(finfo, base):
            print(finfo.name, "is blocked: Hard link to", finfo.linkname)
        elif finfo.islnk() and _badlink(finfo, base):
            print(finfo.name, "is blocked: Symlink to", finfo.linkname)
        else:
            yield finfo


if __name__ == "__main__":
    import doctest

    doctest.testmod()

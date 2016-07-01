"""
Script to set-up the basic directory structure that eqcorrscan likes to use. \
Note that you **do not** have to use this script, or this directory structure \
however it affords a quick-start option for projects, and enforces multiple \
directories for multiple projects.
"""


def _setup_print(msg):
    """
    Print setup specific messages.

    :type msg: str
    """
    prefix = 'EQcorrscan setup:'
    print(' '.join([prefix, msg]))


def create_directory_structure(top_dir):
    """
    Function to create an empty eqcorrscan directory structure.

    :type top_dir: str
    :param top_dir: Top directory to either look into, or create.
    """
    import os

    if not os.path.isdir(top_dir):
        msg = ' '.join([top_dir, 'does not exist: making it.'])
        _setup_print(msg)
        os.makedirs(top_dir)

    for directory in ['parameters', 'templates', 'detections', 'working',
                      'plot']:
        os.makedirs(os.path.join(top_dir, directory))
        msg = ' '.join(['Made', directory, 'directory'])
        _setup_print(msg)

    _setup_print(msg)


def default_matched_filter(top_dir):
    """
    Function to put the default, basic script into the working directory of \
    the project.

    :type top_dir: str
    :param top_dir: Path to eqcorrscan directory structure, must have a \
        working directory, parameter directory and detection directory.
    """
    import os
    import shutil

    if not _is_eqcorrscan_dir(top_dir):
        raise IOError(top_dir + ' is not an eqcorrscan project')
    # Copy script from eqcorrscan install
    default_script = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                  'eqcorrscan_base.py')
    shutil.copyfile(default_script, os.path.join(top_dir, 'working',
                                                 'run_matched_filter.py'))


def _is_eqcorrscan_dir(directory):
    """
    Check to see if the setup is as it should be
    """
    import os
    for required_dir in ['parameters', 'templates', 'detections', 'working']:
        if not os.path.join(directory, required_dir):
            return False
    return True


if __name__ == '__main__':
    import sys
    from eqcorrscan.utils.gui import parameter_setup
    if not len(sys.argv) == 2:
        msg = ('Usage: needs the full or relative path to start an ' +
               'eqcorrscan project within.')
        raise IOError(msg)
    top_dir = sys.argv[1]
    create_directory_structure(top_dir)
    default_matched_filter(top_dir)
    parameter_setup.run()

---
name: Version Checklist
about: Track standard to-do for a release
title: ''
labels: ''
assignees: ''

---

**Its time for a new version of EQcorrscan, first are there any open issues or pull-requests that need to be finished?***


**To-do before release**
- [ ] Check docs;
- [ ] Make sure changes are represented in `CHANGES.md`;
- [ ] Bump version number;
- [ ] Merge develop into master;
- [ ] Release a release candidate on github;
- [ ] Test release candidate on conda-forge: see [pr 4](https://github.com/conda-forge/eqcorrscan-feedstock/pull/4);

If release candidate works on conda-forge then:
- [ ] Generate full release on github;
- [ ] Get DOI of new version from zenodo and update master README.md
- [ ] Release on PyPi (`python setup.py sdist; python -m twine upload dist/*`);
- [ ] Generate PR for eqcorrscan-feedstock from PyPi release;
- [ ] Email google-groups email list.

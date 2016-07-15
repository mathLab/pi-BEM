#!/bin/bash -e

# Settings
REPO_PATH=git@github.com:mathLab/pi-BEM.git
HTML_PATH=../doc/html
COMMIT_USER="Documentation Builder"
COMMIT_EMAIL="luca.heltai@gmail.com"
CHANGESET=$(git rev-parse --verify HEAD)

# Get a clean version of the HTML documentation repo.
rm -rf ${HTML_PATH}
mkdir -p ${HTML_PATH}
git clone -b gh-pages "${REPO_PATH}" --single-branch ${HTML_PATH}

# rm all the files through git to prevent stale files.
cd ${HTML_PATH}
git rm -rf .
cd -

# Generate the HTML documentation.
wget http://www.dealii.org/developer/doxygen/deal.tag
ninja doc

# Create and commit the documentation repo.
cd ${HTML_PATH}
git add .
git config user.name "${COMMIT_USER}"
git config user.email "${COMMIT_EMAIL}"
git commit -m "Automated documentation build for changeset ${CHANGESET}."
git pull --rebase origin gh-pages
git push origin gh-pages
cd -

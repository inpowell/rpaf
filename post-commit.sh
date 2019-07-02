#!/usr/bin/env bash

STASH_NAME=$(cat $GIT_DIR/commit_info/stash)

git stash list | head -n 1 | grep -q $STASH_NAME &&
    git stash pop ||
	true # this is okay

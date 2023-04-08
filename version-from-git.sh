#!/usr/bin/env bash

if [[ "$1" == "-h" ||  "$1" = "--help" ]]; then
	echo "Set the version string from the git currently cloned and checked out."
	exit 0
fi

version="$(git describe --tags|sed -E 's/^v//;s/-([0-9]+)-.*$/.dev\1/')"

if [[ -x "$(which poetry 2> /dev/null)" ]]; then
	exec poetry version "${version}"
else
	exec sed -E -i 's/^(version *= *")[^"]+(")$/\1'"${version}"'\2/' pyproject.toml
fi

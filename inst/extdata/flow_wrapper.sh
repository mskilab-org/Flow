#!/bin/bash -l
[ -z "$CLEANED_SET" ] && exec /usr/bin/env -i CLEANED_SET=1 /bin/bash "$0" "$@"

[ "${CLEANED_SET+set}" = "set" ] && printf "Running in clean environment\n"
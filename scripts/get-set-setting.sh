#!/bin/bash


##
## get or set (if not set) specified setting for a specified setting
##


# script filename
script_name=$(basename "${BASH_SOURCE[0]}")

# check for correct number of arguments
if [ $# -lt 2 ] ; then
	echo -e "\n $script_name ERROR: WRONG NUMBER OF ARGUMENTS SUPPLIED \n" >&2
	echo -e "\n USAGE: $script_name settings_file setting_name [setting_value] \n" >&2
	exit 1
fi

# arguments
settings_txt=$(readlink -f "$1")
setting_name=$2
setting_value_arg=$3


#########################


# check that input exists

if [ ! -s "$settings_txt" ] ; then
	echo -e "\n $script_name ERROR: file $settings_txt does not exist \n" >&2
	exit 1
fi


#########################


# check if setting is already set and output that

setting_value=$(cat "$settings_txt" | grep "^${setting_name}|" | head -1 | cut -d "|" -f 2)
if [ -n "$setting_value" ] ; then
	echo "$setting_value"
	exit
fi


#########################


# determine proper setting value

genome_root=$(cat "$settings_txt" | grep "^GENOME-DIR|" | head -1 | cut -d "|" -f 2)

if [[ "$setting_name" == REF-* ]] ; then
	# reference file settings
	ref_type=${setting_name/REF-/}
	setting_value=$(bash ${BASH_SOURCE%/*}/get-ref.sh "$genome_root" "$ref_type")
elif [ -n "$setting_value_arg" ] ; then
	setting_value=$setting_value_arg
else
	echo -e "\n $script_name ERROR: the value for setting $setting_name must be specified \n" >&2
	exit 1
fi


#########################


# save setting value

if [ -n "$setting_value" ] ; then
	echo "${setting_name}|${setting_value}" >> $settings_txt
else
	echo -e "\n $script_name ERROR: value for $setting_name could not be determined \n" >&2
	exit 1
fi


#########################


# try to get the setting again from the updated settings file

setting_value=$(cat "$settings_txt" | grep "^${setting_name}|" | head -1 | cut -d "|" -f 2)
echo "$setting_value"


#########################



# end

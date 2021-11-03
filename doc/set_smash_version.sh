version=""
# we check if git is available
if [ $(which git) ] 
then
    # we check if a .git directory exists
    if [ -d ${source_dir}/.git ]
    then
        version=$(git describe)
    fi
fi
# if git returned a version, we use it to point to the development tree
if [ $version ]
then
    sed -i "s/smash-devel/smash-devel\/tree\/$version/" \
        ${binary_dir}/user/index.html
else 
    # if git did not return a version, we extract the latest stable version
    # number from the file Changelog.md
    # Note: the case in which the user has downloaded an archive version of the
    # development branch is not taken into account, the link will always be to
    # the last stable version
    tmp=$(grep "## \[SMASH" ${source_dir}/CHANGELOG.md)
    version=${tmp:10:5}
    # if the version is the format x.y.z, z must be a number
    last_digit=${version:4:1}
    if [[ ! $last_digit =~ [[:digit:]] ]]
    then
        # the version is in the format x.y, we remove .z
        version=${version:0:3}
    fi
    version=SMASH-$version
    sed -i "s/smash-devel/smash\/releases\/tag\/$version/" ${binary_dir}/user/index.html
fi

repos="math config predef type_traits static_assert mpl preprocessor assert
    core throw_exception lexical_cast range iterator detail concept_check
    utility numeric_conversion integer array container move smart_ptr"

for repo in $repos
do
    if [ ! -d "$repo" ]; then
        git clone "git@github.com:boostorg/$repo.git"
    fi
done

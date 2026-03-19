process GET_GIT_COMMIT {
    label 'process_single_low'
    label 'git'
 
    input:
    val commit_hash

    output:
    path "git_commit_hash.txt",             emit: commit_hash
    tuple val(task.process),
        val("git"),
        eval('git --version | cut -d" " -f3'), topic: versions
   
    script:
    """
    echo "\nYou can obtain the exact version of the vsearchpipeline used for this run" >> git_commit_hash.txt
    echo "by cloning the repo and running the Git command:\n" >> git_commit_hash.txt
    echo "git checkout ${commit_hash}\n" >> git_commit_hash.txt
    """
}

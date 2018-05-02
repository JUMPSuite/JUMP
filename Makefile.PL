use 5.6.0;
use ExtUtils::MakeMaker;

WriteMakefile(
<<<<<<< HEAD
    EXE_FILES    => [ "bin/jump_sj.pl", "bin/jump_sj_log.pl" ],
    VERSION_FROM => "bin/jump_sj.pl",
=======
    EXE_FILES    => [ "bin/jump_q.pl", "bin/_jump_q.pl" ],
    VERSION_FROM => "bin/jump_q.pl",
>>>>>>> jump-quantification/master
    NAME         => "JUMP",
);

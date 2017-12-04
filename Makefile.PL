use 5.6.0;
use ExtUtils::MakeMaker;

WriteMakefile(
    EXE_FILES    => [ "bin/jump_q.pl" ],
    VERSION_FROM => "bin/jump_q.pl",
    NAME         => "JUMP",
);

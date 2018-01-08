use 5.6.0;
use ExtUtils::MakeMaker;

WriteMakefile(
    EXE_FILES    => [ "bin/jump_f.pl", "bin/_jump_f.pl" ],
    VERSION_FROM => "bin/parallel/jump_f.pl",
    NAME         => "JUMP",
);

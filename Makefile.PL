use 5.6.0;
use ExtUtils::MakeMaker;

WriteMakefile(
    EXE_FILES    => [ "bin/jump_d.pl", "bin/_jump_d.pl" ],
    VERSION_FROM => "bin/_jump_d.pl",
    NAME         => "JUMP",
);

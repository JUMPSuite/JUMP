use 5.6.0;
use ExtUtils::MakeMaker;

WriteMakefile(
    EXE_FILES    => [ "bin/jump_sj.pl", "bin/ReAdW.exe" ],
    VERSION_FROM => "bin/jump_sj.pl",
    NAME         => "JUMP",
);

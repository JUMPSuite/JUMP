use 5.6.0;
use ExtUtils::MakeMaker;

WriteMakefile(
    EXE_FILES    => [ "bin/jump_i.pl"],
    VERSION_FROM => "bin/jump_i.pl",
    NAME         => "JUMPi",
);

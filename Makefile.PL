use 5.6.0;
use ExtUtils::MakeMaker;

WriteMakefile(
    EXE_FILES    => [ "bin/jump_v.pl"],
    VERSION_FROM => "bin/jump_v.pl",
    NAME         => "JUMPv",
);

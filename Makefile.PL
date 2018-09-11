use 5.6.0;
use ExtUtils::MakeMaker;

WriteMakefile(
    EXE_FILES    => [ "bin/jump_aq.pl"],
    VERSION_FROM => "bin/jump_aq.pl",
    NAME         => "JUMPaq",
);

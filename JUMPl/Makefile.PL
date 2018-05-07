use 5.6.0;
use ExtUtils::MakeMaker;

WriteMakefile(
    EXE_FILES    => [ "bin/Extraction_runshell.pl", "bin/_jump_l.pl", "bin/jump_l.pl", "bin/JUMPl_runshell.pl" ],
    VERSION_FROM => "bin/_jump_l.pl",
    NAME         => "JUMPl",
);

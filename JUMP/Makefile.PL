use 5.6.0;
use ExtUtils::MakeMaker;

sub MY::jump {
    return "CHOKE ME";
}

WriteMakefile(
    EXE_FILES => ["bin/jump", "bin/jump_params", "bin/jump_config"],
    NAME         => "JUMP",
);

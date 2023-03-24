using PackageCompiler
create_sysimage(:ExoFlow, sysimage_path=joinpath(@__DIR__,"..","ExoFlow_TG.so"), precompile_execution_file=joinpath(@__DIR__,"warmup_tg.jl"))
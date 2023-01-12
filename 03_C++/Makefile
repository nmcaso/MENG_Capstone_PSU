cpp_files = ConstructDAS.cpp DAS.cpp ImageArea.cpp dataset.cpp
cudafiles = CUDADAS.cu
args_debg = -g -G
args_prod = -O3
lib_incld = -L "C:/Program Files/MATLAB/R2022b/extern/lib/win64/microsoft/" -llibmat -llibmx
compile   = nvcc.exe

all: debug program

debug: $(cpp_files) $(cudafiles)
	$(compile) $(cpp_files) $(cudafiles) -o DAS_dbg.exe $(args_debg) $(lib_incld)

program: $(cpp_files) $(cudafiles)
	$(compile) $(cpp_files) $(cudafiles) -o DAS.exe $(args_prod) $(lib_incld)

run:
	./DAS_dbg.exe

clean:
	rm DAS.exe DAS.lib DAS.pdb DAS_dbg.exe DAS_dbg.lib DAS_dbg.pdb DAS.exp DAS_dbg.exp
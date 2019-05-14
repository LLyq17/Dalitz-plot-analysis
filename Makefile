#--------------------------------

CMAKE_COMMAND = /home/fb/.local/lib/python3.5/site-packages/cmake/data/bin/cmake
CXX=/usr/local/cuda/bin/nvcc
CXXFLAGS=--expt-relaxed-constexpr -Xnvlink=--disable-warnings -Xcompiler=-Wno-attributes --generate-code=arch=compute_61,code=sm_61 --expt-relaxed-constexpr -Xnvlink=--disable-warnings -Xcompiler=-Wno-attributes --generate-code=arch=compute_61,code=sm_61 --expt-relaxed-constexpr -Xnvlink=--disable-warnings -Xcompiler=-Wno-attributes --generate-code=arch=compute_61,code=sm_61 -O2 -g -DNDEBUG -Xcompiler=-fPIC -Wno-deprecated-gpu-targets -shared
CXXLIBS=-L/usr/local/cuda/lib64/stubs  -L/usr/local/cuda/lib64 ../../src/PDFs/basic/libPDFBasic.a ../../src/PDFs/combine/libPDFCombine.a ../../src/PDFs/physics/libPDFPhysics.a ../../src/PDFs/libPDFCore.a ../../src/goofit/libFitManager1.a ../../src/goofit/libApplication.a  ../../src/goofit/libFaddeeva.a ../../src/goofit/libFunctorWriter.a ../../src/goofit/libPdfBase.a ../../src/goofit/libFitManager2.a ../../src/goofit/libDataSet.a ../../src/goofit/libVariable.a ../../extern/FeatureDetector/src/x86/libFeatureDetector.a ../../extern/fmt/fmt/libfmt.a -lcudadevrt -lcudart_static -lrt -lpthread -ldl

CCFLAGS=-Wall -Wextra -Wno-unknown-pragmas -Wno-long-long -Wno-attributes -Wno-sign-compare -Wno-unused-parameter -O2 -g -DNDEBUG 
CCLIBS=-L/usr/local/cuda/lib64/stubs  -L/usr/local/cuda/lib64 -Wl,-rpath,/usr/local/cuda/lib64:/home/fb/root/lib ../../src/PDFs/basic/libPDFBasic.a ../../src/PDFs/combine/libPDFCombine.a ../../src/PDFs/physics/libPDFPhysics.a ../../src/PDFs/libPDFCore.a ../../src/goofit/libFitManager1.a ../../src/goofit/libApplication.a /usr/local/cuda/lib64/libcudart.so ../../src/goofit/libFaddeeva.a ../../src/goofit/libFunctorWriter.a ../../src/goofit/libPdfBase.a ../../src/goofit/libFitManager2.a ../../src/goofit/libDataSet.a ../../src/goofit/libVariable.a /home/fb/root/lib/libCore.so /home/fb/root/lib/libRIO.so /home/fb/root/lib/libNet.so /home/fb/root/lib/libHist.so /home/fb/root/lib/libGraf.so /home/fb/root/lib/libGraf3d.so /home/fb/root/lib/libGpad.so /home/fb/root/lib/libTree.so /home/fb/root/lib/libRint.so /home/fb/root/lib/libPostscript.so /home/fb/root/lib/libMatrix.so /home/fb/root/lib/libPhysics.so /home/fb/root/lib/libMathCore.so /home/fb/root/lib/libThread.so /home/fb/root/lib/libMultiProc.so /home/fb/root/lib/libMinuit.so /home/fb/root/lib/libMinuit2.so -m64 ../../extern/FeatureDetector/src/x86/libFeatureDetector.a ../../extern/fmt/libfmt.a -lcudadevrt -lcudart_static -lrt -lpthread -ldl

ROOT_LIBS = $(shell $(ROOTSYS)/bin/root-config --libs) -lTreePlayer -lMinuit -lXMLIO -lMLP -lRIO -lTMVA


LINKOBJ=cmake_device_link.o
all: DalitzFit_Ds3pi

include flags.make
%.cpp.o:  %.cpp flags.make
	@echo  "\033[32m Building CXX object $@ \033[0m"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o $@ -c $<
	
#DalitzFit_Ds3pi.cpp.o:flags.make DalitzFit_Ds3pi.cpp
#	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o DalitzFit_Ds3pi.cpp.o -c DalitzFit_Ds3pi.cpp

#cmake_device_link.o:DalitzFit_Ds3pi.cpp.o dlink.txt
#	@echo "\033[32m Linking CUDA device code cmake_device_link.o \033[0m"
#	$(CMAKE_COMMAND) -E cmake_link_script dlink.txt --verbose=$(VERBOSE)

#cmake_device_link.o:DalitzFit_Ds3pi.cpp.o dlink.txt
#	@echo "\033[32m Linking CUDA device code cmake_device_link.o \033[0m"
#	$(CXX) $(CXXFALGS) -dlink $< -o $@ $(CXXLIBS)

#DalitzFit_Ds3pi:link.txt DalitzFit_Ds3pi.cpp.o 
#	@echo "\033[32m Linking CXX executable $@ \033[0m"
#	$(CMAKE_COMMAND) -E cmake_link_script link.txt --verbose=$(VERBOSE)
#	@echo "\033[32m $@ done \033[0m"




DalitzFit_Ds3pi: DalitzFit_Ds3pi.cpp.o
	@echo "\033[32m Linking CUDA device code cmake_device_link.o \033[0m"
	$(CXX) $(CXXFALGS) -dlink $< -o $(LINKOBJ) $(CXXLIBS)
	@echo "\033[32m Linking CXX executable $@ \033[0m"
	/usr/bin/c++ $(CCFLAGS)	$<  $(LINKOBJ)  -o $@  $(CCLIBS) 
	@echo "\033[32m $@ done \033[0m"

Dalitz_mipwa: Dalitz_mipwa.cpp.o
	@echo "\033[32m Linking CUDA device code cmake_device_link.o \033[0m"
	$(CXX) $(CXXFALGS) -dlink $< -o $(LINKOBJ) $(CXXLIBS)
	@echo "\033[32m Linking CXX executable $@ \033[0m"
	/usr/bin/c++ $(CCFLAGS) $<  $(LINKOBJ)  -o $@  $(CCLIBS)
	@echo "\033[32m $@ done \033[0m"

dalitzplot: dalitzplot.cpp.o
	@echo "\033[32m Linking CUDA device code cmake_device_link.o \033[0m"
	$(CXX) $(CXXFALGS) -dlink $< -o $(LINKOBJ) $(CXXLIBS)
	@echo "\033[32m Linking CXX executable $@ \033[0m"
	/usr/bin/c++ $(CCFLAGS) $<  $(LINKOBJ)  -o $@  $(CCLIBS)
	@echo "\033[32m $@ done \033[0m"

dalitz_2: dalitz_2.cpp.o
	@echo "\033[32m Linking CUDA device code cmake_device_link.o \033[0m"
	$(CXX) $(CXXFALGS) -dlink $< -o $(LINKOBJ) $(CXXLIBS)
	@echo "\033[32m Linking CXX executable $@ \033[0m"
	/usr/bin/c++ $(CCFLAGS) $<  $(LINKOBJ)  -o $@  $(CCLIBS)
	@echo "\033[32m $@ done \033[0m"

dalitz1: dalitz1.cpp.o
	@echo "\033[32m Linking CUDA device code cmake_device_link.o \033[0m"
	$(CXX) $(CXXFALGS) -dlink $< -o $(LINKOBJ) $(CXXLIBS)
	@echo "\033[32m Linking CXX executable $@ \033[0m"
	/usr/bin/c++ $(CCFLAGS) $<  $(LINKOBJ)  -o $@  $(CCLIBS)
	@echo "\033[32m $@ done \033[0m"


DalitzFit_Ds3pi_pwa:DalitzFit_Ds3pi_pwa.cpp.o
	@echo "\033[32m Linking CUDA device code cmake_device_link.o \033[0m"
	$(CXX) $(CXXFALGS) -dlink $< -o $(LINKOBJ) $(CXXLIBS)
	@echo "\033[32m Linking CXX executable $@ \033[0m"
	/usr/bin/c++ $(CCFLAGS) $<  $(LINKOBJ)  -o $@  $(CCLIBS)
	@echo "\033[32m $@ done \033[0m"


objects=DalitzFit_Ds3pi

clean:
	@rm *.o *.swp $(objects)

.PHONY: clean

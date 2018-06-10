################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../AxisSymmetry.cpp \
../FindRoot.cpp \
../InteriorPoint.cpp \
../PrintSolutions.cpp \
../ReadInput.cpp \
../SlopeWallPoint.cpp \
../WallPoint.cpp \
../main.cpp 

OBJS += \
./AxisSymmetry.o \
./FindRoot.o \
./InteriorPoint.o \
./PrintSolutions.o \
./ReadInput.o \
./SlopeWallPoint.o \
./WallPoint.o \
./main.o 

CPP_DEPS += \
./AxisSymmetry.d \
./FindRoot.d \
./InteriorPoint.d \
./PrintSolutions.d \
./ReadInput.d \
./SlopeWallPoint.d \
./WallPoint.d \
./main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



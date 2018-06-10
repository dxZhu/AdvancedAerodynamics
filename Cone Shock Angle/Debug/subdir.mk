################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../InitialCondition.cpp \
../Iteration.cpp \
../PrintSolution.cpp \
../Properties.cpp \
../main.cpp 

OBJS += \
./InitialCondition.o \
./Iteration.o \
./PrintSolution.o \
./Properties.o \
./main.o 

CPP_DEPS += \
./InitialCondition.d \
./Iteration.d \
./PrintSolution.d \
./Properties.d \
./main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



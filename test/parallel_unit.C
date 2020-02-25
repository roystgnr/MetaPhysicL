#include <metaphysicl/metaphysicl_config.h>

#if defined(METAPHYSICL_HAVE_TIMPI) && defined(METAPHYSICL_HAVE_TIMPI_LIB)

#include <metaphysicl/parallel_dynamicsparsenumberarray.h>
#include <metaphysicl/parallel_dualnumber.h>

#include <timpi/communicator.h>
#include <timpi/parallel_implementation.h>
#include <timpi/timpi_init.h>

#include <vector>
#include <cstdlib>
#include <unistd.h>
#include <cstdio>
#include <iostream>
#include <thread>

#define METAPHYSICL_UNIT_ASSERT(expr)                                                              \
  if (!(expr))                                                                                     \
  metaphysicl_error()

#define METAPHYSICL_UNIT_FP_ASSERT(test, true_value, tolerance)                                    \
  if (!(test < true_value + tolerance && test > true_value - tolerance))               \
  metaphysicl_error()

#define TOLERANCE 1e-12

using namespace TIMPI;
using namespace MetaPhysicL;

Communicator * TestCommWorld;

void
testContainerAllGather()
{
  typedef DualNumber<double, DynamicSparseNumberArray<double, unsigned int>> DualReal;

  std::vector<DualReal> vals;
  const unsigned int my_rank = TestCommWorld->rank();

  // Initialize value
  DualReal in = my_rank;
  // Initialize derivative
  in.derivatives().insert(my_rank) = 1.;

  TestCommWorld->allgather(in, vals);

  const std::size_t comm_size = TestCommWorld->size();
  const std::size_t vec_size = vals.size();

  METAPHYSICL_UNIT_ASSERT(comm_size == vec_size);

  for (std::size_t i = 0; i < vec_size; ++i)
  {
    const auto & dn = vals[i];
    METAPHYSICL_UNIT_FP_ASSERT(dn.value(), double(i), TOLERANCE);
    METAPHYSICL_UNIT_FP_ASSERT(dn.derivatives()[i], double(1), TOLERANCE);
  }
}

int
main(int argc, const char * const * argv)
{
  TIMPI::TIMPIInit init(argc, argv);
  TestCommWorld = &init.comm();

  // volatile int i = 0;
  // char hostname[256];
  // gethostname(hostname, sizeof(hostname));
  // printf("PID %d on %s ready for attach\n", getpid(), hostname);
  // fflush(stdout);
  // while (0 == i)
  //   sleep(5);

  // std::string command = "/usr/bin/lldb";

  // std::cout << "Starting in debugger using: " << command << std::endl;

  // char hostname[256];
  // gethostname(hostname, sizeof(hostname));

  // std::stringstream command_stream;

  // // This will start XTerm and print out some info first... then run the debugger
  // command_stream << "xterm -e \"echo 'Rank: " << TestCommWorld->rank() << "  Hostname: " << hostname
  //                << "  PID: " << getpid() << "'; echo ''; ";

  // // Figure out how to run the debugger
  // if (command.find("lldb") != std::string::npos || command.find("gdb") != std::string::npos)
  //   command_stream << command << " -p " << getpid();

  // // Finish up the command
  // command_stream << "\""
  //                << " & ";

  // std::string command_string = command_stream.str();
  // std::cout << "Running: " << command_string << std::endl;

  // int ret = std::system(command_string.c_str());

  // // Sleep to allow time for the debugger to attach
  // std::this_thread::sleep_for(std::chrono::seconds(10));

  testContainerAllGather();

  return 0;
}

#endif // METAPHYSICL_HAVE_TIMPI && METAPHYSICL_HAVE_TIMPI_LIB

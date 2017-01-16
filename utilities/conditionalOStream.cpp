

#include <conditionalOStream.h>


ConditionalOStream::ConditionalOStream(std::ostream &stream, const bool active)
  :  output_stream (stream),   active_flag(active)
{}


void ConditionalOStream::set_condition(bool flag)
{
  active_flag = flag;
}


bool ConditionalOStream::is_active() const
{
  return active_flag;
}





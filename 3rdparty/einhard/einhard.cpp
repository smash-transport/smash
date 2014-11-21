/**
 * Copyright 2010 Matthias Bach <marix@marix.org>
 * Copyright 2014 Matthias Kretz <kretz@kde.org>
 *
 * This file is part of Einhard.
 *
 * Einhard is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Einhard is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Einhard.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <einhard.hpp>
#include <stdexcept>

namespace einhard
{
char const VERSION[] = "0.4";
#ifndef EINHARD_NO_THREAD_LOCAL
namespace
{
thread_local std::ostringstream t_out;
}  // unnamed namespace
#endif

template <> const char *colorForLogLevel<TRACE>() noexcept
{
	return DBlue_t_::ANSI();
}
template <> const char *colorForLogLevel<DEBUG>() noexcept
{
	return Blue_t_::ANSI();
}
template <> const char *colorForLogLevel<INFO>() noexcept
{
	return DGreen_t_::ANSI();
}
template <> const char *colorForLogLevel<WARN>() noexcept
{
	return Orange_t_::ANSI();
}
template <> const char *colorForLogLevel<ERROR>() noexcept
{
	return DRed_t_::ANSI();
}
template <> const char *colorForLogLevel<FATAL>() noexcept
{
	return Red_t_::ANSI();
}

template <> char const *getLogLevelString<ALL>() noexcept
{
	return "  ALL";
}
template <> char const *getLogLevelString<TRACE>() noexcept
{
	return "TRACE";
}
template <> char const *getLogLevelString<DEBUG>() noexcept
{
	return "DEBUG";
}
template <> char const *getLogLevelString<INFO>() noexcept
{
	return " INFO";
}
template <> char const *getLogLevelString<WARN>() noexcept
{
	return " WARN";
}
template <> char const *getLogLevelString<ERROR>() noexcept
{
	return "ERROR";
}
template <> char const *getLogLevelString<FATAL>() noexcept
{
	return "FATAL";
}
template <> char const *getLogLevelString<OFF>() noexcept
{
	return "  OFF";
}

const char *getLogLevelString( LogLevel level )
{
	switch( level )
	{
	case TRACE:
		return getLogLevelString<TRACE>();
	case DEBUG:
		return getLogLevelString<DEBUG>();
	case INFO:
		return getLogLevelString<INFO>();
	case WARN:
		return getLogLevelString<WARN>();
	case ERROR:
		return getLogLevelString<ERROR>();
	case FATAL:
		return getLogLevelString<FATAL>();
	case ALL:
		return getLogLevelString<ALL>();
	case OFF:
		return getLogLevelString<OFF>();
	}
	// FIXME: better this would throw...
	return "";
}

LogLevel getLogLevel( const std::string &level) {
  if (level == "ALL") {
    return einhard::ALL;
  } else if (level == "TRACE") {
    return einhard::TRACE;
  } else if (level == "DEBUG") {
    return einhard::DEBUG;
  } else if (level == "INFO") {
    return einhard::INFO;
  } else if (level == "WARN") {
    return einhard::WARN;
  } else if (level == "ERROR") {
    return einhard::ERROR;
  } else if (level == "FATAL") {
    return einhard::FATAL;
  } else if (level == "OFF") {
    return einhard::OFF;
  } else {
    throw std::invalid_argument("invalid logging level " + level +
                                ". Accepted values are ALL, TRACE, DEBUG, "
                                "INFO, WARN, ERROR, FATAL, and OFF.");
  }
}

template <LogLevel VERBOSITY> void UnconditionalOutput::doInit( const char *areaName )
{
#ifdef EINHARD_NO_THREAD_LOCAL
	out = &realOut;
#else
	out = &t_out;
#endif
	if( colorize )
	{
		// set color according to log level
		*out << colorForLogLevel<VERBOSITY>();
	}

	// Figure out current time
	time_t rawtime;
	tm *timeinfo;
	time( &rawtime );
	timeinfo = localtime( &rawtime );

	// output it
	const auto oldFill = out->fill();
	out->fill( '0' );
	*out << '[';
	*out << std::setw( 2 ) << timeinfo->tm_hour;
	*out << '\'';
	*out << std::setw( 2 ) << timeinfo->tm_min;
	*out << '\'';
	*out << std::setw( 2 ) << timeinfo->tm_sec;
	*out << ']';
	out->fill( oldFill );
	// TODO would be good to have this at least .01 seconds
	// for non-console output pure timestamp would probably be better

	// output the log level and logging area of the message
	*out << ' ' << getLogLevelString<VERBOSITY>();
	if( areaName && areaName[0] != '\0' )
	{
		*out << ' ' << areaName;
	}
	*out << ": ";

	indent = out->str().size();
	if( colorize )
	{
		// The bytes from the ANSI color code don't appear on screen and thus must be subtracted from the indent
		// value. We still make an error if the areaName contains multi-byte (utf-8) characters. At
		// this point that's just not supported.
		indent -= sizeof( "\33[00;30m" ) - 1;  // sizeof includes the trailing \0
		*out << NoColor_t_::ANSI();
	}
}

template void UnconditionalOutput::doInit<TRACE>( const char * );
template void UnconditionalOutput::doInit<DEBUG>( const char * );
template void UnconditionalOutput::doInit<INFO>( const char * );
template void UnconditionalOutput::doInit<WARN>( const char * );
template void UnconditionalOutput::doInit<ERROR>( const char * );
template void UnconditionalOutput::doInit<FATAL>( const char * );

void UnconditionalOutput::checkColorReset()
{
	if( resetColor )
	{
		*out << NoColor_t_::ANSI();
		resetColor = false;
	}
}

void UnconditionalOutput::doCleanup(std::FILE *outfile) noexcept
{
	*out << '\n';
	std::string s = out->str();
	out->str( std::string() );
	std::string s2;
	s2.reserve( s.size() + 18 * 2 );
	std::size_t pos = 0;
	std::size_t start = 0;
	while( ( pos = s.find( '\n', start ) + 1 ) <
	       s.size() - 1 )  // the size - 1 catches double \n\n at the end of the string, reducing it to a single \n
	{
		s2.append( s, start, pos - start );
		s2.append( indent, ' ' );
		start = pos;
	}
	s2.append( s, start, pos - start );
	std::fwrite( s2.c_str(), s2.size(), 1, outfile );
	std::fflush( outfile );  // FIXME: don't want to flush too often, what's the right logic here?
}
}  // namespace einhard

// vim: ts=4 sw=4 tw=100 noet

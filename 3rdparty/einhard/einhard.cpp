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

namespace einhard
{
char const VERSION[] = "0.4";
thread_local std::ostringstream t_out;

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

template <LogLevel VERBOSITY> void OutputFormatter::doInit( const char *areaName )
{
	out = &t_out;
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
	*out << '[';
	*out << std::setfill( '0' ) << std::setw( 2 ) << timeinfo->tm_hour;
	*out << ':';
	*out << std::setfill( '0' ) << std::setw( 2 ) << timeinfo->tm_min;
	*out << ':';
	*out << std::setfill( '0' ) << std::setw( 2 ) << timeinfo->tm_sec;
	*out << ']';
	// TODO would be good to have this at least .01 seconds
	// for non-console output pure timestamp would probably be better

	// output the log level and logging area of the message
	*out << ' ' << getLogLevelString<VERBOSITY>();
	if( areaName )
	{
		*out << ' ' << areaName;
	}
	*out << ": ";

	if( colorize )
	{
		*out << NoColor_t_::ANSI();
	}
}

template void OutputFormatter::doInit<TRACE>( const char * );
template void OutputFormatter::doInit<DEBUG>( const char * );
template void OutputFormatter::doInit<INFO>( const char * );
template void OutputFormatter::doInit<WARN>( const char * );
template void OutputFormatter::doInit<ERROR>( const char * );
template void OutputFormatter::doInit<FATAL>( const char * );

void OutputFormatter::checkColorReset()
{
	if( resetColor )
	{
		*out << NoColor_t_::ANSI();
		resetColor = false;
	}
}

void OutputFormatter::doCleanup() noexcept
{
	*out << '\n';
	const std::string s = out->str();
	out->str( std::string() );
	std::fwrite( s.c_str(), s.size(), 1, stdout );
}
}  // namespace einhard

// vim: ts=4 sw=4 tw=100 noet

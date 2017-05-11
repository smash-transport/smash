/**
 * @file
 *
 * This is the main include file for Einhard.
 *
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

/** \page Einhard Logging Library
 *
 * \section intro_sec Introduction
 *
 * Einhard aims at being a lightweight logging library with the following features:
 *  - Severity filtering
 *  - Automatic output colorization
 *  - Output timestamping
 *  - Low performance overhead
 *
 * To use Einhard all you need to do is create an einhard::Logger object and push messages into its output
 * streams.
 *
 * \code
 * using namespace einhard;
 * Logger logger( INFO );
 * logger.trace() << "Trace message"; // will not be printed
 * logger.error() << "Error message"; // will be printed
 * \endcode
 *
 * To reduce the performance overhad you can specify NDEBUG during compile. This will disable trace
 * and debug messages in a way that should allow the compilers dead code elimination to remove
 * everything that is only pushed into such a stream.
 *
 * \section install_sec Installation
 *
 * Einhard is build using cmake. You can install it using the usual cmake triplet:
 * \code
 * cmake
 * make
 * sudo make install
 * \endcode
 * If you want to build a static lib or install to a custom path you can use the usual cmake
 * utilities to adjust the configuration.
 */

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <cstring>
#include <sstream>
#include <bitset>

#include "stacktrace.h"

// This C header is sadly required to check whether writing to a terminal or a file
#include <cstdio>

#include <unistd.h>  // for isatty

#ifdef __GNUC__
#define EINHARD_ALWAYS_INLINE_ __attribute__((always_inline))
#else
#define EINHARD_ALWAYS_INLINE_
#endif

// Error on MacOS:
//    “thread-local storage is unsupported for the current target”
// To enable a workaround the EINHARD_NO_THREAD_LOCAL macro must be defined.
#ifdef __APPLE__
#define EINHARD_NO_THREAD_LOCAL 1
#endif

/**
 * This namespace contains all objects required for logging using Einhard.
 */
namespace einhard
{
	/**
	 * Version string of the Einhard library
	 */
	extern char const VERSION[];

	/**
	 * Specification of the message severity.
	 *
	 * In output each level automatically includes the higher levels.
	 */
	enum LogLevel
	{
		ALL,   /**< Log all message */
		TRACE, /**< The lowes severity for messages describing the program flow */
		DEBUG, /**< Debug messages */
		INFO,  /**< Messages of informational nature, expected processing time e.g. */
		WARN,  /**< Warning messages */
		ERROR, /**< Non-fatal errors */
		FATAL, /**< Messages that indicate terminal application failure */
		OFF    /**< If selected no messages will be output */
	};

	/**
	 * Retrieve a human readable representation of the given log level value.
	 *
	 * The overload can optimize better because it can determine the LogLevel at compile time.
	 */
	template <LogLevel> const char *getLogLevelString() noexcept;
	/**
	 * Overload of the above function for situations where the LogLevel \p level is only determined at run time.
	 */
	const char *getLogLevelString( LogLevel level );

	/**
	 * Compares the string \p level against the strings for LogLevel and returns the one it matches.
	 *
	 * \param level A string, which is a textual representation of one of the
	 *              LogLevel enumerators.
	 * \return The enumerator that matches the input string.
	 * \throws std::invalid_argument if the string does not match any enumerator.
	 */
	LogLevel getLogLevel( const std::string &level );

	template <LogLevel> const char *colorForLogLevel() noexcept;

	/**
	 * A stream modifier that allows to colorize the log output.
	 */
	template <typename Parent> class Color
	{
	public:
		/// The default color modifier only affects the next object in the stream.
		EINHARD_ALWAYS_INLINE_ Color() noexcept : reset( true )
		{
		}
		/// With the ~ operator the color modifier affects the rest of the stream (or until
		/// another color object is used).
		EINHARD_ALWAYS_INLINE_ Color<Parent> operator~() const noexcept
		{
			return {false};
		}
		EINHARD_ALWAYS_INLINE_ char const *ansiCode() const noexcept
		{
			return Parent::ANSI();
		}
		EINHARD_ALWAYS_INLINE_ bool resetColor() const noexcept
		{
			return reset;
		}

	private:
		EINHARD_ALWAYS_INLINE_ Color( bool r ) noexcept : reset( r )
		{
		}
		bool reset;
	};
#define _COLOR( name, code )                                                                                           \
	struct name##_t_                                                                                               \
	{                                                                                                              \
		static char const *ANSI() noexcept                                                                     \
		{                                                                                                      \
			return "\33[" code "m";                                                                        \
		}                                                                                                      \
	};                                                                                                             \
	typedef Color<name##_t_> name

	_COLOR(DGray,   "01;30");
	_COLOR(Black,   "00;30");
	_COLOR(Red,     "01;31");
	_COLOR(DRed,    "00;31");
	_COLOR(Green,   "01;32");
	_COLOR(DGreen,  "00;32");
	_COLOR(Yellow,  "01;33");
	_COLOR(Orange,  "00;33");
	_COLOR(Blue,    "01;34");
	_COLOR(DBlue,   "00;34");
	_COLOR(Magenta, "01;35");
	_COLOR(DMagenta,"00;35");
	_COLOR(Cyan,    "01;36");
	_COLOR(DCyan,   "00;36");
	_COLOR(White,   "01;37");
	_COLOR(Gray,    "00;37");
	_COLOR(NoColor, "0"    );
#undef _COLOR

	/**
	 * A minimal class that implements the output stream operator to do nothing. This completely
	 * eliminates the output stream statements from the resulting binary.
	 */
	struct DummyOutputFormatter
	{
		template <typename T> EINHARD_ALWAYS_INLINE_ DummyOutputFormatter &operator<<( const T & ) noexcept
		{
			return *this;
		}
		EINHARD_ALWAYS_INLINE_ DummyOutputFormatter &operator<<(
		    std::ostream &( * )( std::ostream & ) ) noexcept
		{
			return *this;
		}
	};

	class UnconditionalOutput
	{
	private:
		// Pointer to the thread_local stringstream (if enabled)
		std::ostringstream *out;
#ifdef EINHARD_NO_THREAD_LOCAL
		// without thread_local we simply use a local stringstream object
		std::ostringstream realOut;
#endif
		// The number of chars required for aligning
		unsigned char indent;
		// Whether to colorize the output
		const bool colorize;
		// Whether the color needs to be reset with the next operator<<
		bool resetColor = false;

	public:
		template <LogLevel VERBOSITY>
		EINHARD_ALWAYS_INLINE_ UnconditionalOutput( bool colorize_, const char *areaName,
							    std::integral_constant<LogLevel, VERBOSITY> )
		    : colorize( colorize_ )
		{
			doInit<VERBOSITY>( areaName );
		}

		template <typename T> UnconditionalOutput &operator<<( const Color<T> &col )
		{
			if( colorize )
			{
				*out << col.ansiCode();
				resetColor = col.resetColor();
			}
			return *this;
		}

		EINHARD_ALWAYS_INLINE_ UnconditionalOutput &operator<<( std::ostream &( *manip )( std::ostream & ) )
		{
			*out << manip;
			return *this;
		}

		template <typename T> EINHARD_ALWAYS_INLINE_ UnconditionalOutput &operator<<( const T &msg )
		{
			*out << msg;
			checkColorReset();
			return *this;
		}

		void doCleanup(std::FILE *outfile = stdout) noexcept;

	protected:
		EINHARD_ALWAYS_INLINE_ UnconditionalOutput( bool colorize_ ) : colorize( colorize_ )
		{
		}
		template <LogLevel VERBOSITY> void doInit( const char *areaName );
		void checkColorReset();
	};
	/**
	 * A wrapper for the output stream taking care proper formatting and colorization of the output.
	 */
	class OutputFormatter : public UnconditionalOutput
	{
		private:
			// the output file (normally used to switch from stdout to stderr)
			std::FILE *const outfile;
			// Whether output is enabled
			const bool enabled;

		public:
			OutputFormatter( const OutputFormatter & ) = delete;
			OutputFormatter( OutputFormatter && ) = default;

			template <LogLevel VERBOSITY>
			EINHARD_ALWAYS_INLINE_ OutputFormatter( bool enabled_, bool const colorize_,
								const char *areaName,
								std::integral_constant<LogLevel, VERBOSITY>,
								std::FILE *outfile_ = stdout )
			    : UnconditionalOutput( colorize_ ), enabled( enabled_ ), outfile( outfile_ )
			{
				if( enabled )
				{
					doInit<VERBOSITY>( areaName );
				}
			}

			template <typename T> OutputFormatter &operator<<( const Color<T> &col )
			{
				if( enabled )
				{
					UnconditionalOutput::operator<<(col);
				}
				return *this;
			}

			EINHARD_ALWAYS_INLINE_ OutputFormatter &operator<<( std::ostream &( *manip )( std::ostream & ) )
			{
				if( enabled )
				{
					UnconditionalOutput::operator<<(manip);
				}
				return *this;
			}

			template <typename T> EINHARD_ALWAYS_INLINE_ OutputFormatter &operator<<( const T &msg )
			{
				if( enabled )
				{
					UnconditionalOutput::operator<<(msg);
				}
				return *this;
			}

			EINHARD_ALWAYS_INLINE_
			~OutputFormatter( )
			{
				if( enabled )
				{
					doCleanup( outfile );
				}
			}
	};

	/**
     * A Logger object can be used to output messages to stdout.
     *
     * The Logger object is created with a certain verbosity. All messages of a more verbose
     * LogLevel will be filtered out. The way the class is build this can happen at compile
     * time up to the level restriction given by the template parameter.
     *
     * The class can automatically detect non-tty output and will not colorize output in that case.
     */ 
	template<LogLevel MAX = ALL> class Logger
	{
		private:
			char areaName[32 - sizeof( LogLevel ) - sizeof( bool )] = {'\0'};
			LogLevel verbosity;
			bool colorize_stdout;
			bool colorize_stderr;

		public:
			/**
			 * Create a new Logger object.
			 *
			 * The object will automatically colorize output on ttys and not colorize output
			 * on non ttys.
			 */
			Logger( const LogLevel verbosity = WARN ) : verbosity( verbosity )
			{
				// use some, sadly not c++-ways to figure out whether we are writing ot a terminal
				// only colorize when we are writing ot a terminal
				colorize_stdout = isatty( fileno( stdout ) );
				colorize_stderr = isatty( fileno( stderr ) );
			};
			/**
			 * Create a new Logger object explicitly selecting whether to colorize the output or not.
			 *
			 * Be aware that if output colorization is selected output will even be colorized if
			 * output is to a non tty.
			 */
			Logger( const LogLevel verbosity, const bool colorize )
			    : verbosity( verbosity ), colorize_stdout( colorize ), colorize_stderr( colorize ) {};

			/**
			 * Set an area name. This will be printed after the LogLevel to identify the
			 * place in the code where the output is coming from. This can be used to
			 * identify the different Logger objects in the log output.
			 *
			 * \param name A string. Only the first 30, or so, characters will be used. The rest
			 *             will not be displayed. You can reset the name with an empty string.
			 * \warning Passing a nullptr is not allowed!
			 */
			void setAreaName( const char *name )
			{
				assert( name );
				std::strncpy( &areaName[0], name, sizeof( areaName ) - 1 );
				areaName[sizeof( areaName ) - 1] = '\0';
			}
			EINHARD_ALWAYS_INLINE_
			void setAreaName( const std::string &name )
			{
				setAreaName(name.c_str());
			}

			/** Access to the trace message stream. */
#ifdef NDEBUG
			DummyOutputFormatter trace() const noexcept
			{
				return DummyOutputFormatter();
			}

			template <typename... Ts> void trace( Ts &&... ) const noexcept
			{
			}
#else
			OutputFormatter trace() const
			{
				return {isEnabled<TRACE>(), colorize_stdout, areaName,
					std::integral_constant<LogLevel, TRACE>()};
			}

			template <typename... Ts> void trace( Ts &&... args ) const noexcept
			{
				if( isEnabled<TRACE>() )
				{
					UnconditionalOutput o{colorize_stdout, areaName,
							      std::integral_constant<LogLevel, TRACE>()};
					auto &&unused = {&( o << args )...};
					o.doCleanup();
				}
			}
#endif
			/** Access to the debug message stream. */
#ifdef NDEBUG
			DummyOutputFormatter debug() const noexcept
			{
				return DummyOutputFormatter();
			}
			template <typename... Ts> void debug( Ts &&... ) const noexcept
			{
			}
#else
			OutputFormatter debug() const
			{
				return {isEnabled<DEBUG>(), colorize_stdout, areaName,
					std::integral_constant<LogLevel, DEBUG>()};
			}
			template <typename... Ts> void debug( Ts &&... args ) const noexcept
			{
				if( isEnabled<DEBUG>() )
				{
					UnconditionalOutput o{colorize_stdout, areaName,
							  std::integral_constant<LogLevel, DEBUG>()};
					auto &&unused = {&( o << args )...};
					o.doCleanup();
				}
			}
#endif
			/** Access to the info message stream. */
			OutputFormatter info() const
			{
				return {isEnabled<INFO>(), colorize_stdout, areaName,
					std::integral_constant<LogLevel, INFO>()};
			}
			template <typename... Ts> void info( Ts &&... args ) const noexcept
			{
				if( isEnabled<INFO>() )
				{
					UnconditionalOutput o{colorize_stdout, areaName,
							  std::integral_constant<LogLevel, INFO>()};
					auto &&unused = {&( o << args )...};
					o.doCleanup();
				}
			}
			/** Access to the warning message stream. */
			OutputFormatter warn() const
			{
				return {isEnabled<WARN>(), colorize_stderr, areaName,
					std::integral_constant<LogLevel, WARN>(), stderr};
			}
			template <typename... Ts> void warn( Ts &&... args ) const noexcept
			{
				if( isEnabled<WARN>() )
				{
					UnconditionalOutput o{colorize_stderr, areaName,
							  std::integral_constant<LogLevel, WARN>()};
					auto &&unused = {&( o << args )...};
					o.doCleanup(stderr);
				}
			}
			/** Access to the error message stream. */
			OutputFormatter error() const
			{
				return {isEnabled<ERROR>(), colorize_stderr, areaName,
					std::integral_constant<LogLevel, ERROR>(), stderr};
			}
			template <typename... Ts> void error( Ts &&... args ) const noexcept
			{
				if( isEnabled<ERROR>() )
				{
					UnconditionalOutput o{colorize_stderr, areaName,
							  std::integral_constant<LogLevel, ERROR>()};
					auto &&unused = {&( o << args )...};
					o.doCleanup(stderr);
				}
			}
			/** Access to the fatal message stream. */
			OutputFormatter fatal() const
			{
				return {isEnabled<FATAL>(), colorize_stderr, areaName,
					std::integral_constant<LogLevel, FATAL>(), stderr};
			}
			/** Will print stacktrace. */
			template <typename... Ts> void fatal( Ts &&... args ) const noexcept
			{
				if( isEnabled<FATAL>() )
				{
					UnconditionalOutput o{colorize_stderr, areaName,
							  std::integral_constant<LogLevel, FATAL>()};
					auto &&unused = {&( o << args )...};
					o.doCleanup(stderr);
					print_stacktrace(stderr, 63, 2);
				}
			}

			template <LogLevel LEVEL> bool isEnabled() const noexcept
			{
#ifdef NDEBUG
				if( LEVEL == DEBUG || LEVEL == TRACE ) {
					return false;
				}
#endif
				return ( MAX <= LEVEL && verbosity <= LEVEL );
			}

			/** Modify the verbosity of the Logger.
 			 *
			 * Be aware that the verbosity can not be increased over the level given by the template
			 * parameter
			 */
			inline void setVerbosity( LogLevel verbosity ) noexcept
			{
				this->verbosity = verbosity;
			}
			/** Retrieve the current log level.
			 *
			 * If you want to check whether messages of a certain LogLevel will be output the
			 * method isEnabled() should be
			 * preffered.
			 */
			inline LogLevel getVerbosity() const noexcept
			{
				return this->verbosity;
			}
			/**
			 * Retrieve a human readable representation of the current log level
			 */
			inline char const * getVerbosityString( ) const
			{
				return getLogLevelString(verbosity);
			}
			/**
			 * Select whether the output stream should be colorized.
			 */
			void setColorize( bool colorize ) noexcept
			{
				this->colorize_stdout = colorize;
				this->colorize_stderr = colorize;
			}
			/**
			 * Check whether the output stream is colorized.
			 */
			bool getColorize() const noexcept
			{
				return this->colorize_stdout || this->colorize_stderr;
			}
	};
}

// vim: ts=4 sw=4 tw=100 noet

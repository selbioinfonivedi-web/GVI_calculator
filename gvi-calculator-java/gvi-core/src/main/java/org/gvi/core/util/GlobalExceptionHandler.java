package org.gvi.core.util;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Crash containment (Section 2): install this once at application start
 * (CLI main() and, later, the JavaFX Application.start()) so that any
 * exception which escapes every other catch block still gets logged with a
 * full stack trace and reported to the user, instead of the JVM dumping to
 * stderr and dying silently or with a raw trace the user can't act on.
 */
public final class GlobalExceptionHandler {

    private static final Logger log = LoggerFactory.getLogger(GlobalExceptionHandler.class);

    private GlobalExceptionHandler() {
    }

    public interface CrashListener {
        void onUncaughtCrash(Thread thread, Throwable throwable);
    }

    public static void install(CrashListener listener) {
        Thread.setDefaultUncaughtExceptionHandler((thread, throwable) -> {
            log.error("Uncaught exception on thread '{}' -- this is a bug, please report it with the log file attached",
                    thread.getName(), throwable);
            if (listener != null) {
                try {
                    listener.onUncaughtCrash(thread, throwable);
                } catch (Throwable t) {
                    log.error("Crash listener itself threw", t);
                }
            }
        });
    }

    public static void install() {
        install(null);
    }
}

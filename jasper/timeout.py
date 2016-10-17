from __future__ import with_statement
import signal
from contextlib import contextmanager

'''
http://stackoverflow.com/questions/366682/how-to-limit-execution-time-of-a-function-call-in-python
'''


class TimeoutException(Exception):
    pass


@contextmanager
def time_limit(seconds):
    def signal_handler(signum, frame):
        raise TimeoutException("Timed out!")
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)

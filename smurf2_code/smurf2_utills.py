
import time
import logging


def time_it(method):
    # logging.info('Start {}'.format(method.__name__))
    def timed(*args, **kw):
        logging.info('Entered %r' % method.__name__)
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        t = te - ts
        if t > 60:
            logging.info('Done %s: %3.2f minutes\n' % (method.__name__, t/60))
        else:
            logging.info('Done %s: %3.2f sec\n' % (method.__name__, t))
        return result
    return timed


def define_logger(level=logging.INFO, file_name=None):
    log_level = level
    FORMAT = '%(asctime)-15s %(levelname)s %(funcName)s() [%(lineno)s] - %(message)s'
    logging.basicConfig(format=FORMAT, level=log_level, filename=file_name)

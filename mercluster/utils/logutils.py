import logging
import os

def getLogger(name, logfile):
	"""
	Create a logger that writes warnings and above to the console,
	and debug and above to a file

	Args:
		name: unique name for logger
		logfile: file to use for logging
	Returns:
		logger object
	"""
	logger = logging.getLogger(name)
	logger.setLevel(logging.DEBUG)
	ch = logging.StreamHandler()
	ch.setLevel(logging.WARNING)
	formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
	ch.setFormatter(formatter)

	fh = logging.FileHandler(logfile)
	fh.setLevel(logging.DEBUG)
	formatter = logging.Formatter(
		'%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	fh.setFormatter(formatter)
	logger.addHandler(ch)
	logger.addHandler(fh)

	return logger

def closeLoggerHandlers(logger):
	"""
	Closes out the files and handlers in a logger object
	Args:
		logger: the logger to close handlers from
	"""
	handlers = logger.handlers[:]
	for handler in handlers:
		handler.close()
		logger.removeHandler(handler)

def restart_log(logfile: str):
	"""
	Deletes a log file and replaces it with an empty file

	Args:
		logfile: path to log file
	"""
	os.remove(logfile)


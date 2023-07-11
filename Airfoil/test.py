from model.util import MyLogger


if __name__ == '__main__':
    my_logger1=MyLogger(log_filename='test1.log')
    my_logger2=MyLogger(log_filename='test2.log')
    my_logger1.logger.info('1')
    my_logger2.logger.info('2')

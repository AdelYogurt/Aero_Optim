classdef Logger < handle
    % create a logger object
    %
    properties
        log_filename;
        loglevel;
        formatter;
    end
    methods
        function self=Logger(log_filename, loglevel)
            if nargin < 2
                loglevel='INFO';
            end
            % create logger
            self.log_filename=log_filename;
            self.loglevel=loglevel;

            % default create a file log
            self.info('initialize logger');
        end

        function set.loglevel(self,loglevel)
            % set log level
            %
            switch loglevel
                case 'DEBUG'
                    digtial_loglevel=1;
                case 'INFO'
                    digtial_loglevel=2;
                case 'WARNING'
                    digtial_loglevel=3;
                case 'ERROR'
                    digtial_loglevel=4;
                case 'CRITICAL'
                    digtial_loglevel=5;
                otherwise
            end
            self.loglevel=digtial_loglevel;
        end

        function debug(self,message)
            if self.loglevel > 0
                file_handler=fopen(self.log_filename,'a');
                fprintf(file_handler,['[DEBUG]',datestr(now),' %s\n'],message);
                fclose(file_handler);
                clear('file_handler');
            end
        end

        function info(self,message)
            if self.loglevel > 1
                file_handler=fopen(self.log_filename,'a');
                fprintf(file_handler,['[INFO]',datestr(now),' %s\n'],message);
                fclose(file_handler);
                clear('file_handler');
            end
        end

        function warning(self,message)
            if self.loglevel > 2
                file_handler=fopen(self.log_filename,'a');
                fprintf(file_handler,['[WARNING]',datestr(now),' %s\n'],message);
                fclose(file_handler);
                clear('file_handler');
            end
        end

        function error(self,message)
            if self.loglevel > 3
                file_handler=fopen(self.log_filename,'a');
                fprintf(file_handler,['[ERROR]',datestr(now),' %s\n'],message);
                fclose(file_handler);
                clear('file_handler');
            end
        end

        function critial(self,message)
            if self.loglevel > 4
                file_handler=fopen(self.log_filename,'a');
                fprintf(file_handler,['[CRITIAL]',datestr(now),' %s\n'],message);
                fclose(file_handler);
                clear('file_handler');
            end
        end
    end
end
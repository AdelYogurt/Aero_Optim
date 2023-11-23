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
            log_filedir=fileparts(log_filename);
            if ~isfolder(log_filedir)
                mkdir(log_filedir)
            end
            self.info('initialize logger');
        end

        function set.loglevel(self,loglevel)
            % set log level
            %
            switch loglevel
                case {'DEBUG',0}
                    digtial_loglevel=0;
                case {'INFO',1}
                    digtial_loglevel=1;
                case {'WARNING',2}
                    digtial_loglevel=2;
                case {'ERROR',3}
                    digtial_loglevel=3;
                case {'CRITICAL',4}
                    digtial_loglevel=4;
                otherwise
            end
            self.loglevel=digtial_loglevel;
        end

        function debug(self,message)
            if self.loglevel < 1
                file_handler=fopen(self.log_filename,'a');
                fprintf(file_handler,['[DEBUG]',datestr(now),' ']);fprintf(file_handler,'%s\n',message);
                fclose(file_handler);
                clear('file_handler');
            end
        end

        function info(self,message)
            if self.loglevel < 2
                file_handler=fopen(self.log_filename,'a');
                fprintf(file_handler,['[INFO]',datestr(now),' ']);fprintf(file_handler,'%s\n',message);
                fclose(file_handler);
                clear('file_handler');
            end
        end

        function warning(self,message)
            if self.loglevel < 3
                file_handler=fopen(self.log_filename,'a');
                fprintf(file_handler,['[WARNING]',datestr(now),' ']);fprintf(file_handler,'%s\n',message);
                fclose(file_handler);
                clear('file_handler');
            end
        end

        function error(self,message)
            if self.loglevel < 4
                file_handler=fopen(self.log_filename,'a');
                fprintf(file_handler,['[ERROR]',datestr(now),' ']);fprintf(file_handler,'%s\n',message);
                fclose(file_handler);
                clear('file_handler');
            end
        end

        function critial(self,message)
            if self.loglevel < 5
                file_handler=fopen(self.log_filename,'a');
                fprintf(file_handler,['[CRITIAL]',datestr(now),' ']);fprintf(file_handler,'%s\n',message);
                fclose(file_handler);
                clear('file_handler');
            end
        end
    end
end
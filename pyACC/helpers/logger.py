from datetime import datetime
from typing import Any
from colorama import Fore, Style
import time 

class Logger:
    def __init__(self, alias: str):
        self.alias = alias

    def __call__(self, message):
        kind = "INFO"
        msg = f"{datetime.fromtimestamp(time.time())} :: {Style.BRIGHT + Fore.MAGENTA + self.alias + Style.RESET_ALL} :: {Fore.GREEN + kind + Style.RESET_ALL} :: {message}"

        print(msg)

    def warning(self, message):
        kind = "WARNING"
        msg = f"{datetime.fromtimestamp(time.time())} :: {Style.BRIGHT + Fore.MAGENTA + self.alias + Style.RESET_ALL} :: {Fore.YELLOW + kind + Style.RESET_ALL} :: {message}"

        print(msg)

    def error(self, message, ErrorType=None):
        kind = "ERROR"
        msg = f"{datetime.fromtimestamp(time.time())} :: {Style.BRIGHT + Fore.MAGENTA + self.alias + Style.RESET_ALL} :: {Fore.RED + kind + Style.RESET_ALL} :: {message}"

        if ErrorType is None:
            print(msg)
        else:
            raise(ErrorType(msg))

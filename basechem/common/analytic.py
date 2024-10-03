import logging

# Events that can be logged
PAGE_VIEW = "page_view"
LOGIN = "login"
TASK = "task"


class Analytic:
    """
    This class exists to capture analytics in AWS cloudwatch. Creating an Analytic object
    sends a relevant logging statement to the "analytics" logger, as configured in the django
    settings
    """

    def __init__(self, event, project_code, username, **kwargs):
        self.event = event
        self.project_code = project_code
        self.username = username
        self.other_params = kwargs
        self.log()

    def log(self):
        """
        Log this analytic to the analytics logger
        """
        logging_statement = f"[{self.event}] [{self.project_code}] [{self.username}]"

        if self.other_params:
            for param, value in self.other_params.items():
                logging_statement += f" [{param}={value}]"

        analytics_logger = logging.getLogger("analytics")
        analytics_logger.info(logging_statement)

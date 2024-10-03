from django.contrib import messages
from django.contrib.auth.models import AbstractUser
from django.db import models

from basechem.common.admin_emailer import send_easter_egg_email

# Keeps track of how many points users need to earn to find an easter egg
REQUIRED_EASTER_EGG_POINTS = {"Hiking": 3, "ESP": 1, "Torsion": 3, "Assay MMPs": 5}


class BasechemUser(AbstractUser):
    # A dictionary mapping the name of an egg hunt (ex. "Hiking") to the number
    # of points the user has scored towards that hunt
    easter_egg_points = models.JSONField(default=dict, blank=True, null=True)

    def __str__(self):
        return self.get_full_name()

    def get_full_name(self):
        """Return the first and last names, minus any spaces at end or beginning."""
        return f"{self.first_name} {self.last_name}".strip()

    def get_short_name(self):
        """Return the first name, if it is not blank; otherwise, return last name."""
        return self.first_name if self.first_name else self.last_name

    def get_email_name(self):
        """
        Return the user's email without the domain name. Return the entire email address
        if there is no '@domain_name'
        """
        if "@" in self.email:
            return self.email.partition("@")[0]
        else:
            return self.email

    def add_easter_egg_point(self, http_request, egg_name):
        """
        Adds a point to the user's `easter_egg_points` for `egg_name`. If the user has reached the required
        number of points to find the easter egg, sends an email to the user and MnI to award
        a prize/brownie points
        :param http_request: the http request from the view that called this function (to show a message in the browser)
        :param egg_name: a string, the name of the easter egg whose points should be incremented
        """
        # Add an easter egg point
        if not self.easter_egg_points.get(egg_name):
            self.easter_egg_points[egg_name] = 0
        self.easter_egg_points[egg_name] += 1
        self.save()

        # Award the easter egg if the user has enough points
        if self.easter_egg_points[egg_name] == REQUIRED_EASTER_EGG_POINTS[egg_name]:
            messages.info(
                http_request,
                f"Woooo! You've found the {egg_name} Easter egg! Check your email for a surprise :)",
            )
            send_easter_egg_email(self, egg_name)

from django.conf import settings
from django.core.mail import send_mail
from django.template import loader


def send_easter_egg_email(user, egg_name):
    """
    Sends an email to MnI and the user who found an Easter egg
    :param user: the user who found the Easter egg
    :param egg_name: the name of the Easter egg that was found
    """
    html_message = loader.render_to_string(
        "main/components/easter_egg_email.html", {"user": user, "egg_name": egg_name}
    )
    send_mail(
        "Basechem Easter Egg Found!",
        f"{user.first_name} has found the '{egg_name}' Easter egg!",
        settings.DEFAULT_FROM_EMAIL,
        [settings.ADMIN_EMAIL, user.email],
        html_message=html_message,
    )

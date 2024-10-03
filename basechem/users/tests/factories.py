import factory

from basechem.users.models import BasechemUser


class BasechemUserFactory(factory.DjangoModelFactory):
    first_name = factory.Faker("first_name")
    last_name = factory.Faker("last_name")
    username = factory.Faker("user_name")
    email = factory.Faker("email")
    password = factory.PostGenerationMethodCall("set_password", "testpassword")

    class Meta:
        model = BasechemUser

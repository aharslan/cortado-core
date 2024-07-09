import multiprocessing
import threading


class Singleton:
    def __init__(self, decorated):
        self._decorated = decorated
        self._lock = threading.Lock()
        self._instance = None

    def instance(self):
        if self._instance is None:
            with self._lock:
                # double check because other thread could have generated the instance in the meantime
                if self._instance is None:
                    self._instance = self._decorated()

        return self._instance

    def __call__(self):
        raise TypeError("Singletons must be accessed through `instance()`.")

    def __instancecheck__(self, inst):
        return isinstance(inst, self._decorated)


@Singleton
class PoolFactory:
    def __init__(self):
        self.pool = multiprocessing.Pool(8)

    def get_pool(self):
        return self.pool

    def restart_pool(self):
        self.pool.terminate()
        self.pool.join()

        self.pool = multiprocessing.Pool()

#define HASHLIN_INIT(name)                      \
  do {                                          \
    name = malloc(sizeof *name);                \
    PANIC_MEM(name, stderr);                    \
    tommy_hashlin_init(name);                   \
  } while (0)

#define HASHLIN_DONE(name, free_func)                                   \
  do {                                                                  \
    tommy_hashlin_foreach(name, (tommy_foreach_func*)free_func);        \
    tommy_hashlin_done(name);                                           \
    free(name);                                                         \
  } while (0)

#define ARRAY_INIT(name)                        \
  do {                                          \
    name = malloc(sizeof *name);                \
    PANIC_MEM(name, stderr);                    \
    tommy_array_init(name);                     \
  } while (0)

#define ARRAY_DONE(name, free_func)                     \
  do {                                                  \
    for (int i = 0; i < tommy_array_size(name); ++i) {  \
      free_func(tommy_array_get(name, i));              \
    }                                                   \
                                                        \
    tommy_array_done(name);                             \
    free(name);                                         \
  } while (0)

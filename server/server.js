const express = require('express');
const bodyParser = require('body-parser');
const mongoose = require('mongoose');
const cookieParser = require('cookie-parser');
const config = require('./config/environment');

const mainRouter = require('./routes/main');
const apiRouter = require('./routes/api');
const authRouter = require('./routes/auth');

const app = express();
const PORT = config.port;

// Connect to MongoDB
mongoose.connect(config.mongo.uri)
.then(() => console.log('Connected to MongoDB'))
.catch(err => console.error('Error connecting to MongoDB:', err));

// Middleware
app.use(bodyParser.json());
app.use(cookieParser());
app.use(express.static('public'));


// Use router
app.use('/', mainRouter);
app.use('/api', apiRouter);
app.use('/auth', authRouter);


app.listen(PORT, () => {
  console.log(`Server running on http://localhost:${PORT}`);
});
